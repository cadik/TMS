from hdr2ldr import HDR2LDR_finetune
import tensorflow as tf
import pprint
import os

flags = tf.app.flags
flags.DEFINE_integer("epoch", 100, "Number of epoch [100]")
flags.DEFINE_integer("batch_size", 32, "The size of batch images [32]")
flags.DEFINE_integer("im_height", 512, "The height of image to use [64]")
flags.DEFINE_integer("im_width", 512, "The size of label to produce [128]")
flags.DEFINE_string("dataset_dir", "dataSet_1", "Name of dataset directory [dataSet_1]")
flags.DEFINE_string("checkpoint_dir", "checkpoint", "Name of checkpoint directory [checkpoint]")
flags.DEFINE_string("sample_dir", "samples", "Name of sample directory for testing [samples]")
flags.DEFINE_string("summary_dir", 'summary', "Name of summary directory [summary]")
flags.DEFINE_string("weights_dir", "newWeights.pkl", "Name of ldr2hdr weights [newWeights.pkl]")
flags.DEFINE_boolean("is_train", False, "True for training, False for testing [True]")
FLAGS = flags.FLAGS

pp = pprint.PrettyPrinter()

def main(_):
    #pp.pprint(flags.FLAGS.__flags)

    if not os.path.exists(FLAGS.checkpoint_dir):
        os.makedirs(FLAGS.checkpoint_dir)
    if not os.path.exists(FLAGS.sample_dir):
        os.makedirs(FLAGS.sample_dir)
    if not os.path.exists(FLAGS.summary_dir):
        os.makedirs(FLAGS.summary_dir)
    with tf.Session() as sess:
        hdr2ldr_finetune = HDR2LDR_finetune(sess,
                        checkpoint_dir=FLAGS.checkpoint_dir,
                        sample_dir=FLAGS.sample_dir,
                        weights_dir=FLAGS.weights_dir,
                        is_train=FLAGS.is_train,
                        batch_size=FLAGS.batch_size,
                        im_height=FLAGS.im_height,
                        im_width=FLAGS.im_width,
                        learning_rate=0.1)

        hdr2ldr_finetune.run(FLAGS)

if __name__=='__main__':
    tf.app.run()