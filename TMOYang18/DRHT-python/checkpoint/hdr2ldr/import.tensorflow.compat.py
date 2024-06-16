import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()

# Load the TensorFlow model
with tf.Session() as sess:
    # Restore the saved model
    saver = tf.train.import_meta_graph("hdr2ldr.model-400.meta")
    saver.restore(sess, "hdr2ldr.model-400")

    # Get the graph
    graph = tf.get_default_graph()

    # Get the names of output tensors
    output_tensor_names = [tensor.name for tensor in tf.get_default_graph().as_graph_def().node if "" in tensor.name]

    # Print the output tensor names
    print("Output tensor names:", output_tensor_names)