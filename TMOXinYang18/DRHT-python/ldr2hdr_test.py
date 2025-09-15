import os, sys
import tensorflow as tf
import tensorlayer as tl
import numpy as np
import ldr2hdr, img_io



FLAGS = tf.flags.FLAGS
tf.flags.DEFINE_integer("width", "512", "Reconstruction image width")
tf.flags.DEFINE_integer("height", "512", "Reconstruction image height")
tf.flags.DEFINE_string("im_dir", "./input", "Path to image directory or an individual image")
tf.flags.DEFINE_string("out_dir", "./output", "Path to output directory")
tf.flags.DEFINE_string("params", "./checkpoint/ldr2hdr.npz", "Path to trained CNN weights")


sx = int(np.maximum(32, np.round(FLAGS.width/32.0)*32))
sy = int(np.maximum(32, np.round(FLAGS.height/32.0)*32))


frames = [FLAGS.im_dir]
if os.path.isdir(FLAGS.im_dir):
    frames = [os.path.join(FLAGS.im_dir, name)
              for name in sorted(os.listdir(FLAGS.im_dir))
              if os.path.isfile(os.path.join(FLAGS.im_dir, name))]

x = tf.placeholder(tf.float32, shape=[1, sy, sx, 3])
net = ldr2hdr.model(x)
y = ldr2hdr.get_final(net, x)
sess = tf.InteractiveSession()

load_params = tl.files.load_npz(name=FLAGS.params)
tl.files.assign_params(sess, load_params, net)

if not os.path.exists(FLAGS.out_dir):
    os.makedirs(FLAGS.out_dir)
k = 0
for i in range(len(frames)):

    
    x_buffer = img_io.readLDR(frames[i], (sy,sx), True, 1.0)
    y_predict = sess.run([y], feed_dict={x: x_buffer})
    y_gamma = np.power(np.maximum(y_predict, 0.0), 0.5)
    k += 1;
    outname = './output' + frames[i][7:-4] + '.exr'
    img_io.writeEXR(y_predict, outname)

sess.close()

