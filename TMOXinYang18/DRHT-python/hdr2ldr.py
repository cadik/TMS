import img_io
import time
import os
import tensorflow as tf
import numpy as np
import cv2
import pickle

class HDR2LDR_finetune(object):
    def __init__(self, sess=None, checkpoint_dir=None,sample_dir=None, weights_dir=None, is_train=False, batch_size=32, im_height=64, im_width=128, learning_rate=0.1):
        self.im_height = im_height
        self.im_width = im_width
        self.batch_size = batch_size
        self.sess = sess
        self.checkpoint_dir = checkpoint_dir
        self.sample_dir = sample_dir
        self.weights_dir = None#weights_dir
        self.is_train = is_train
        self.counter = 0
        self.last_improvement = 0
        self.best_validation = 100.
        self.learning_rate = learning_rate
        self.build_model()
    
    def build_model(self):
        self.images = tf.placeholder(tf.float32, [None, self.im_height, self.im_width, 3], name='images')
        self.labels = tf.placeholder(tf.float32, [None, self.im_height, self.im_width, 3], name='labels')

        with tf.variable_scope('encoder_weights'):
            self.enWeights = {'enW1': tf.get_variable('w1_xavier',[1,1,3,64], initializer=tf.contrib.layers.xavier_initializer(), trainable=True),
                              'enW2': tf.get_variable('w2_xavier',[7,7,64,64], initializer=tf.contrib.layers.xavier_initializer(), trainable=True),
                              'enW3': tf.get_variable('w3_xavier',[5,5,64,128],initializer=tf.contrib.layers.xavier_initializer(), trainable=True),
                              'enW4': tf.get_variable('w4_xavier',[3,3,128,256],initializer=tf.contrib.layers.xavier_initializer(), trainable=True),
                              'enW5': tf.get_variable('w5_xavier',[3,3,256,256],initializer=tf.contrib.layers.xavier_initializer(), trainable=True),
                              'enW6': tf.get_variable('w6_xavier',[3,3,256,512],initializer=tf.contrib.layers.xavier_initializer(), trainable=True)
                            }
        with tf.variable_scope('encoder_biases'):
            self.enBiases = {'enB1': tf.Variable(tf.zeros([64]), trainable=True, name='en_B1'),
                             'enB2': tf.Variable(tf.zeros([64]), trainable=True, name='en_B2'),
                             'enB3': tf.Variable(tf.zeros([128]), trainable=True, name='en_B3'),
                             'enB4': tf.Variable(tf.zeros([256]), trainable=True, name='en_B4'),
                             'enB5': tf.Variable(tf.zeros([256]), trainable=True, name='en_B5'),
                             'enB6': tf.Variable(tf.zeros([512]), trainable=True, name='en_B6'),
                            }


        with tf.variable_scope('decoder_weights'):
            self.deWeights = {'deW1': tf.get_variable('w1_xavier',[3,3,256,512],initializer=tf.contrib.layers.xavier_initializer(), trainable=True),
                              'deW2': tf.get_variable('w2_xavier',[3,3,256,256],initializer=tf.contrib.layers.xavier_initializer(), trainable=True),
                              'deW3': tf.get_variable('w3_xavier',[3,3,128,256],initializer=tf.contrib.layers.xavier_initializer(), trainable=True),
                              'deW4': tf.get_variable('w4_xavier',[5,5,64,128],initializer=tf.contrib.layers.xavier_initializer(), trainable=True),
                              'deW5': tf.get_variable('w5_xavier',[7,7,64,64],initializer=tf.contrib.layers.xavier_initializer(), trainable=True),
                              'deW6': tf.get_variable('w6_xavier',[1,1,64,3],initializer=tf.contrib.layers.xavier_initializer(), trainable=True)
                            }
        with tf.variable_scope('decoder_biases'):
            self.deBiases = {'deB1': tf.Variable(tf.zeros([256]), trainable=True, name='de_B1'),
                             'deB2': tf.Variable(tf.zeros([256]), trainable=True, name='de_B2'),
                             'deB3': tf.Variable(tf.zeros([128]), trainable=True, name='de_B3'),
                             'deB4': tf.Variable(tf.zeros([64]), trainable=True, name='de_B4'),
                             'deB5': tf.Variable(tf.zeros([64]), trainable=True, name='de_B5'),
                             'deB6': tf.Variable(tf.zeros([3]), trainable=True, name='de_B6')
                            }
        self.pred = self.forward() # prediction
        self.saver = tf.train.Saver()

    def run(self, config):
        print('Testing...')
        if self.load(self.checkpoint_dir):
            print(' [*] Load SUCCESS')
        inNames = os.listdir('./output')
        # read on batch for testing, otherwise the performance drops
        for i in xrange(self.batch_size):
            if i == 0:
                inp = img_io.readEXR('./output/' + inNames[i])
                inp = np.expand_dims(inp, 0)
            else:
                inp1 = img_io.readEXR('./output/' + inNames[i])
                inp1 = np.expand_dims(inp1, 0)
                inp = np.concatenate([inp, inp1], 0)
        print inp.shape
        images = self.sess.run(self.pred, feed_dict={self.images:inp})
        for i in xrange(self.batch_size):

            res = np.reshape(images[i,...], [self.im_height, self.im_width,3])
            res = np.minimum(np.maximum(res, 0.), 1.)
            resName = './samples/'+ inNames[i][:-4]+'.jpg'
            cv2.imwrite(resName, res* 255.)


    def forward(self):
        self.expInput = tf.pow(self.images, 1./2.2) / 30.
        self.logInput = tf.log(self.images + 1. / 255.)
        conv1 = tf.nn.conv2d(self.logInput, self.enWeights['enW1'], strides=[1,1,1,1], padding='SAME') + self.enBiases['enB1']
        conv1_bn = tf.contrib.layers.batch_norm(conv1)
        self.conv1_bn_act = tf.nn.elu(conv1_bn)
        conv2 = tf.nn.conv2d(self.conv1_bn_act, self.enWeights['enW2'], strides=[1,2,2,1], padding='SAME') + self.enBiases['enB2']
        conv2_bn = tf.contrib.layers.batch_norm(conv2)
        self.conv2_bn_act = tf.nn.elu(conv2_bn)
        conv3 = tf.nn.conv2d(self.conv2_bn_act, self.enWeights['enW3'], strides=[1,2,2,1], padding='SAME') + self.enBiases['enB3']
        conv3_bn = tf.contrib.layers.batch_norm(conv3)
        self.conv3_bn_act = tf.nn.elu(conv3_bn)
        conv4 = tf.nn.conv2d(self.conv3_bn_act, self.enWeights['enW4'], strides=[1,2,2,1], padding='SAME') + self.enBiases['enB4']
        conv4_bn = tf.contrib.layers.batch_norm(conv4)
        self.conv4_bn_act = tf.nn.elu(conv4_bn)
        conv5 = tf.nn.conv2d(self.conv4_bn_act, self.enWeights['enW5'], strides=[1,2,2,1], padding='SAME') + self.enBiases['enB5']
        conv5_bn = tf.contrib.layers.batch_norm(conv5)
        self.conv5_bn_act = tf.nn.elu(conv5_bn)
        feat = tf.nn.conv2d(self.conv5_bn_act, self.enWeights['enW6'], strides=[1,2,2,1], padding='SAME') + self.enBiases['enB6']
        feat_bn = tf.contrib.layers.batch_norm(feat)
        self.feat_bn_act = tf.nn.elu(feat_bn)
        deconv1 = tf.nn.conv2d_transpose(self.feat_bn_act, self.deWeights['deW1'], output_shape=tf.stack([tf.shape(self.feat_bn_act)[0],tf.shape(self.feat_bn_act)[1]*2,tf.shape(self.feat_bn_act)[2]*2,256]), strides=[1,2,2,1], padding="SAME") + self.deBiases['deB1']
        deconv1_bn = tf.contrib.layers.batch_norm(deconv1)
        self.deconv1_bn_act = tf.nn.elu(deconv1_bn)
        skip_1 = self.deconv1_bn_act + self.conv5_bn_act  + self.deconv1_bn_act * self.conv5_bn_act
        deconv2 = tf.nn.conv2d_transpose(skip_1, self.deWeights['deW2'], output_shape=tf.stack([tf.shape(self.feat_bn_act)[0],tf.shape(skip_1)[1]*2,tf.shape(skip_1)[2]*2,256]), strides=[1,2,2,1], padding="SAME") + self.deBiases['deB2']
        deconv2_bn = tf.contrib.layers.batch_norm(deconv2)
        self.deconv2_bn_act = tf.nn.elu(deconv2_bn)
        skip_2 = self.deconv2_bn_act + self.conv4_bn_act  + self.deconv2_bn_act * self.conv4_bn_act
        deconv3 = tf.nn.conv2d_transpose(skip_2, self.deWeights['deW3'], output_shape=tf.stack([tf.shape(self.feat_bn_act)[0],tf.shape(skip_2)[1]*2,tf.shape(skip_2)[2]*2,128]), strides=[1,2,2,1], padding="SAME") + self.deBiases['deB3']
        deconv3_bn = tf.contrib.layers.batch_norm(deconv3)
        self.deconv3_bn_act = tf.nn.elu(deconv3_bn)
        skip_3 = self.deconv3_bn_act + self.conv3_bn_act  + self.deconv3_bn_act * self.conv3_bn_act
        deconv4 = tf.nn.conv2d_transpose(skip_3, self.deWeights['deW4'], output_shape=tf.stack([tf.shape(self.feat_bn_act)[0],tf.shape(skip_3)[1]*2,tf.shape(skip_3)[2]*2,64]), strides=[1,2,2,1], padding="SAME") + self.deBiases['deB4']
        deconv4_bn = tf.contrib.layers.batch_norm(deconv4)
        self.deconv4_bn_act = tf.nn.elu(deconv4_bn)
        skip_4 = self.deconv4_bn_act + self.conv2_bn_act  + self.deconv4_bn_act * self.conv2_bn_act
        deconv5 = tf.nn.conv2d_transpose(skip_4, self.deWeights['deW5'], output_shape=tf.stack([tf.shape(self.feat_bn_act)[0],tf.shape(skip_4)[1]*2,tf.shape(skip_4)[2]*2,64]), strides=[1,2,2,1], padding="SAME") + self.deBiases['deB5']
        deconv5_bn = tf.contrib.layers.batch_norm(deconv5)
        self.deconv5_bn_act = tf.nn.elu(deconv5_bn)
        skip_5 = self.deconv5_bn_act + self.conv1_bn_act  + self.deconv5_bn_act * self.conv1_bn_act
        resImg = tf.nn.conv2d(skip_5, self.deWeights['deW6'], strides=[1,1,1,1], padding="SAME") + self.deBiases['deB6']
        resImg_bn = tf.contrib.layers.batch_norm(resImg)
        self.resImg_bn_act = tf.nn.elu(resImg_bn + self.expInput)
        return self.resImg_bn_act

    def load(self, checkpoint_dir):
        print(' [*] Reading checkpoints...')
        model_dir = '%s' % ('hdr2ldr')
        checkpoint_dir = os.path.join(checkpoint_dir, model_dir)

        ckpt = tf.train.get_checkpoint_state(checkpoint_dir)
        if ckpt and ckpt.model_checkpoint_path:
            ckpt_name = os.path.basename(ckpt.model_checkpoint_path)
            self.saver.restore(self.sess, os.path.join(checkpoint_dir, ckpt_name))
            return True
        else:
            return False
