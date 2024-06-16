import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()

def convert_to_pb(saved_model_dir, output_pb_file):
    # Load the saved model
    with tf.Session(graph=tf.Graph()) as sess:
        # Restore the saved model
        meta_graph_def = tf.saved_model.loader.load(sess, [tf.saved_model.tag_constants.SERVING], saved_model_dir)

        # Convert variables to constants
        output_graph_def = tf.graph_util.convert_variables_to_constants(
            sess,
            meta_graph_def.graph_def,
            ["output_tensor_names"]  # Add the names of the output tensors you want to freeze
        )

        # Serialize the graph def protocol buffer and write it to a binary file
        with tf.gfile.GFile(output_pb_file, "wb") as f:
            f.write(output_graph_def.SerializeToString())

# Example usage
saved_model_dir = "path/to/saved_model"
output_pb_file = "output_graph.pb"
convert_to_pb(saved_model_dir, output_pb_file)


# net_in = tl.layers.InputLayer(x_in, name='input_layer') hdr2ldr.model-400.index