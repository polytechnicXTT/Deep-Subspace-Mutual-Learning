import tensorflow as tf
import numpy as np
import math
from libs.activations import lrelu
from libs.utils import corrupt


def data_load(data_name):
    data_dir = "./data/COAD/"
    # data_name = "BREAST_Gene_Expression.txt"
    data_path = data_dir + data_name
    f = open(data_path)  # 读取文件
    lines = f.readlines()
    # print(len(lines))

    for i in range(len(lines)):
        lines[i] = lines[i].strip().split("\t")[1:]
        # print(lines[i])

    gene_expression = np.zeros([len(lines) - 1, len(lines[1])])

    for i in range(len(lines)):
        if i == 0:
            continue
        for j in range(len(lines[1])):
            gene_expression[i - 1][j] = float(lines[i][j])

    gene_expression = np.matrix(gene_expression)
    # print(gene_expression.shape)
    gene_expression = np.transpose(gene_expression)
    gene_expression = np.array(gene_expression)

    return gene_expression


# %%
def autoencoder(gene_expression,
                n_filters=[15],
                filter_sizes=[5],
                corruption=False):
    input_shape = [None, gene_expression.shape[1]]
    # input to the network
    x = tf.placeholder(tf.float32, input_shape, name='x')

    # %%
    # ensure 2-d is converted to square tensor.
    if len(x.get_shape()) == 2:
        x_tensor = tf.reshape(
            x, [-1, 1, gene_expression.shape[1], 1])
    elif len(x.get_shape()) == 4:
        x_tensor = x
    else:
        raise ValueError('Unsupported input dimensions')
    current_input = x_tensor

    # %%
    # Optionally apply denoising autoencoder
    if corruption:
        current_input = corrupt(current_input)

    # %%
    # Build the encoder
    encoder = []
    shapes = []
    for layer_i, n_output in enumerate(n_filters[:]):
        n_input = current_input.get_shape().as_list()[3]
        shapes.append(current_input.get_shape().as_list())
        W = tf.Variable(
            tf.random_uniform([
                1,
                filter_sizes[layer_i],
                n_input, n_output],
                -1.0 / math.sqrt(n_input),
                1.0 / math.sqrt(n_input)))
        b = tf.Variable(tf.zeros([n_output]))
        encoder.append(W)
        output = lrelu(
            tf.add(tf.nn.conv2d(
                current_input, W, strides=[1, 1, 2, 1], padding='SAME'), b))
        current_input = output
    enc_w=W
    enc_b=b
    # %%
    # store the latent representation
    z = current_input
    encoder.reverse()
    shapes.reverse()

    # %%
    # Build the decoder using the same weights
    for layer_i, shape in enumerate(shapes):
        W = encoder[layer_i]
        b = tf.Variable(tf.zeros([W.get_shape().as_list()[2]]))
        output = lrelu(tf.add(
            tf.nn.conv2d_transpose(
                current_input, W,
                tf.stack([tf.shape(x)[0], shape[1], shape[2], shape[3]]),
                strides=[1, 1, 2, 1], padding='SAME'), b))
        current_input = output

    dec_w=W
    dec_b=b
    # %%
    # now have the reconstruction through the network
    y = current_input
    # cost function measures pixel-wise difference
    cost = tf.reduce_sum(tf.square(y - x_tensor))

    # %%
    return {'x': x, 'z': z, 'y': y, 'cost': cost,'W_e':enc_w,'W_d':dec_w,'b_e':enc_b,'b_d':dec_b}


def test_mnist(gene_expression):
    """Test the convolutional autoencder using MNIST."""
    pretrain_model_dir = "./pre_train_model/COAD/partial_pretrain/0.5/"

    ae = autoencoder(gene_expression)

    learning_rate = 0.01
    optimizer = tf.train.AdamOptimizer(learning_rate).minimize(ae['cost'])

    # %%
    # We create a session to use the graph
    saver = tf.train.Saver()
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    sess = tf.Session(config=config)
    sess.run(tf.global_variables_initializer())
    # saver.restore(sess, "D:/model1/BREAST/Gene_Expression2410.ckpt")

    n_epochs = 2000
    for epoch_i in range(n_epochs):

        batch_xs = gene_expression[:]
        train = np.array(batch_xs)
        sess.run(optimizer, feed_dict={ae['x']: train})
        print(epoch_i, sess.run(ae['cost'], feed_dict={ae['x']: train}))
        if epoch_i % 10 == 0:
            saver.save(sess, pretrain_model_dir + "Gene_Expression%d.ckpt" % epoch_i)         # 保存pretrain_model

        if epoch_i > 500:
            k = 0.002
            learning_rate = learning_rate / (1 + k * epoch_i)
    W_e = sess.run(ae['W_e'], feed_dict={ae['x']: train})
    print(W_e)
    return{"W_encoder": W_e}


# %%
if __name__ == '__main__':
    test_mnist(data_load("COLON_Gene_0.5.txt"))
