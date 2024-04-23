import tensorflow as tf
import tensorlayer as tl
import numpy as np

def model(x):

    x_in = tf.scalar_mul(255.0, x)
    net_in = tl.layers.InputLayer(x_in, name='input_layer')
    conv_layers, skip_layers = encoder(net_in)

    network = tl.layers.Conv2dLayer(conv_layers,
                    act = tf.identity,
                    shape = [3, 3, 512, 512],
                    strides = [1, 1, 1, 1],
                    padding='SAME',
                    name ='encoder/h6/conv')
    network = tl.layers.BatchNormLayer(network, is_train=False, name='encoder/h6/batch_norm')
    network.outputs = tf.nn.relu(network.outputs, name='encoder/h6/relu')

    network = decoder(network, skip_layers)

    return network

def get_final(network, x_in):
    sb, sy, sx, sf = x_in.get_shape().as_list()
    y_predict = network.outputs
    thr = 0.05
    alpha = tf.reduce_max(x_in, reduction_indices=[3])
    alpha = tf.minimum(1.0, tf.maximum(0.0, alpha-1.0+thr)/thr)
    alpha = tf.reshape(alpha, [-1, sy, sx, 1])
    alpha = tf.tile(alpha, [1,1,1,3])

    x_lin = tf.pow(x_in, 2.0)
    y_predict = tf.exp(y_predict)-1.0/255.0
    y_final = (1-alpha)*x_lin + alpha*y_predict
    
    return y_final


def encoder(input_layer):
    VGG_MEAN = [103.939, 116.779, 123.68]
    red, green, blue = tf.split(input_layer.outputs, 3, 3)
    bgr = tf.concat([ blue - VGG_MEAN[0], green - VGG_MEAN[1], red - VGG_MEAN[2] ], axis=3)
    network = tl.layers.InputLayer(bgr, name='encoder/input_layer_bgr')
    network     = conv_layer(network, [ 3, 64], 'encoder/h1/conv_1')
    beforepool1 = conv_layer(network, [64, 64], 'encoder/h1/conv_2')
    network     = pool_layer(beforepool1, 'encoder/h1/pool')
    network     = conv_layer(network, [64, 128], 'encoder/h2/conv_1')
    beforepool2 = conv_layer(network, [128, 128], 'encoder/h2/conv_2')
    network     = pool_layer(beforepool2, 'encoder/h2/pool')
    network     = conv_layer(network, [128, 256], 'encoder/h3/conv_1')
    network     = conv_layer(network, [256, 256], 'encoder/h3/conv_2')
    beforepool3 = conv_layer(network, [256, 256], 'encoder/h3/conv_3')
    network     = pool_layer(beforepool3, 'encoder/h3/pool')
    network     = conv_layer(network, [256, 512], 'encoder/h4/conv_1')
    network     = conv_layer(network, [512, 512], 'encoder/h4/conv_2')
    beforepool4 = conv_layer(network, [512, 512], 'encoder/h4/conv_3')
    network     = pool_layer(beforepool4, 'encoder/h4/pool')
    network     = conv_layer(network, [512, 512], 'encoder/h5/conv_1')
    network     = conv_layer(network, [512, 512], 'encoder/h5/conv_2')
    beforepool5 = conv_layer(network, [512, 512], 'encoder/h5/conv_3')
    network     = pool_layer(beforepool5, 'encoder/h5/pool')

    return network, (input_layer, beforepool1, beforepool2, beforepool3, beforepool4, beforepool5)


def decoder(input_layer, skip_layers, batch_size=1, is_training=False):
    sb, sx, sy, sf = input_layer.outputs.get_shape().as_list()
    alpha = 0.0
    network = deconv_layer(input_layer, (batch_size,sx,sy,sf,sf), 'decoder/h1/decon2d', alpha, is_training)
    network = skip_connection_layer(network, skip_layers[5], 'decoder/h2/fuse_skip_connection', is_training)
    network = deconv_layer(network, (batch_size,2*sx,2*sy,sf,sf), 'decoder/h2/decon2d', alpha, is_training)
    network = skip_connection_layer(network, skip_layers[4], 'decoder/h3/fuse_skip_connection', is_training)
    network = deconv_layer(network, (batch_size,4*sx,4*sy,sf,sf/2), 'decoder/h3/decon2d', alpha, is_training)
    network = skip_connection_layer(network, skip_layers[3], 'decoder/h4/fuse_skip_connection', is_training)
    network = deconv_layer(network, (batch_size,8*sx,8*sy,sf/2,sf/4), 'decoder/h4/decon2d', alpha, is_training)
    network = skip_connection_layer(network, skip_layers[2], 'decoder/h5/fuse_skip_connection', is_training)
    network = deconv_layer(network, (batch_size,16*sx,16*sy,sf/4,sf/8), 'decoder/h5/decon2d', alpha, is_training)
    network = skip_connection_layer(network, skip_layers[1], 'decoder/h6/fuse_skip_connection', is_training)
    network = tl.layers.Conv2dLayer(network,
                        act = tf.identity,
                        shape = [1, 1, int(sf/8), 3],
                        strides=[1, 1, 1, 1],
                        padding='SAME',
                        W_init = tf.contrib.layers.xavier_initializer(uniform=False),
                        b_init = tf.constant_initializer(value=0.0),
                        name ='decoder/h7/conv2d')
    network = tl.layers.BatchNormLayer(network, is_train=is_training, name='decoder/h7/batch_norm')
    network.outputs = tf.maximum(alpha*network.outputs, network.outputs, name='decoder/h7/leaky_relu')
    network = skip_connection_layer(network, skip_layers[0], 'decoder/h7/fuse_skip_connection')

    return network



def conv_layer(input_layer, sz, str):
    network = tl.layers.Conv2dLayer(input_layer,
                    act = tf.nn.relu,
                    shape = [3, 3, sz[0], sz[1]],
                    strides = [1, 1, 1, 1],
                    padding = 'SAME',
                    name = str)
    return network

def pool_layer(input_layer, str):
    network = tl.layers.PoolLayer(input_layer,
                    ksize=[1, 2, 2, 1],
                    strides=[1, 2, 2, 1],
                    padding='SAME',
                    pool = tf.nn.max_pool,
                    name = str)
    return network

def skip_connection_layer(input_layer, skip_layer, str, is_training=True):
    _, sx, sy, sf = input_layer.outputs.get_shape().as_list()
    _, sx_, sy_, sf_ = skip_layer.outputs.get_shape().as_list()
    assert (sx_,sy_,sf_) == (sx,sy,sf)
    skip_layer.outputs = tf.log(tf.pow(tf.scalar_mul(1.0/255, skip_layer.outputs), 2.0)+1.0/255.0)
    network = tl.layers.ConcatLayer(layer = [input_layer,skip_layer], concat_dim=3, name ='%s/skip_connection'%str)
    network = tl.layers.Conv2dLayer(network,
                    act = tf.identity,
                    shape = [1, 1, sf+sf_, sf],
                    strides = [1, 1, 1, 1],
                    padding = 'SAME',
                    b_init = tf.constant_initializer(value=0.0),
                    name = str)
    return network

def deconv_layer(input_layer, sz, str, alpha, is_training=True):
    scale = 2
    filter_size = (2 * scale - scale % 2)
    num_in_channels = int(sz[3])
    num_out_channels = int(sz[4])
    network = tl.layers.DeConv2dLayer(input_layer,
                                shape = [filter_size, filter_size, num_out_channels, num_in_channels],
                                output_shape = [sz[0], sz[1]*scale, sz[2]*scale, num_out_channels],
                                strides = [1, scale, scale, 1],
                                padding = 'SAME',
                                act = tf.identity,
                                name = str)
    network = tl.layers.BatchNormLayer(network, is_train=is_training, name='%s/batch_norm_dc'%str)
    network.outputs = tf.maximum(alpha*network.outputs, network.outputs, name='%s/leaky_relu_dc'%str)
    return network