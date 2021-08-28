import tensorflow as tf
import numpy as np
from pre_train import data_load, autoencoder
from load_data import load_data
import pandas as pd
import os

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # or any {'0', '1', '2'}
tf.logging.set_verbosity(tf.logging.ERROR)

data_dir = "./data/KRCCC/"
data_name1 = "KIDNEY_Gene_Expression.txt"
data_name2 = "KIDNEY_Methy_Expression.txt"
data_name3 = "KIDNEY_Mirna_Expression.txt"

pretrain_model_dir = "./pre_train_model/KRCCC/"
train_model_dir = "./train_model/KRCCC/"
xlsx_dir = "./model/KRCCC/"


def test_mnist(gene_expression):
    """Test the convolutional autoencder using MNIST."""
    # load MNIST as before
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
    saver.restore(sess, pretrain_model_dir + "Gene_Expression1990.ckpt")         #  加载预训练模型

    # batch_size = gene_expression.shape[0]
    n_epochs = 21
    for epoch_i in range(n_epochs):

        batch_xs = gene_expression[:]
        train = np.array(batch_xs)
        sess.run(optimizer, feed_dict={ae['x']: train})
        print(epoch_i, sess.run(ae['cost'], feed_dict={ae['x']: train}))
        if epoch_i % 10 == 0:
            saver.save(sess, train_model_dir + "model_test1_%d.ckpt" % epoch_i)          # 保存训练模型

        # if epoch_i > 500:
        #     k = 0.002
        #     learning_rate = learning_rate / (1 + k * epoch_i)
    W_e, W_d, b_e, b_d = sess.run([ae['W_e'], ae['W_d'], ae['b_e'], ae['b_d']], feed_dict={ae['x']: train})
    print(W_e)
    return{"W_encoder": W_e, "W_decoder": W_d, "bias_e": b_e, "bias_d": b_d}


class ConvAE:
    def __init__(self, n_input, kernel_size, n_hidden, reg_const1=1.0, reg_const2=150.0, batch_size=256, denoise=False):
        # n_hidden is a arrary contains the number of neurals on every layer
        self.n_input = n_input
        self.n_hidden = n_hidden
        self.kernel_size = kernel_size
        self.iter = 0
        self.batch_size = batch_size
        net_weights = self._initialize_weights()
        self.weights = net_weights

        # model_1
        self.x1 = tf.placeholder(tf.float32, [batch_size, self.n_input[0][0], self.n_input[0][1], 1])
        self.x2 = tf.placeholder(tf.float32, [batch_size, self.n_input[1][0], self.n_input[1][1], 1])
        self.x3 = tf.placeholder(tf.float32, [batch_size, self.n_input[2][0], self.n_input[2][1], 1])

        self.learning_rate = tf.placeholder(tf.float32, [])

        if denoise is False:
            x_input1 = self.x1
            x_input2 = self.x2
            x_input3 = self.x3
        else:
            x_input1 = tf.add(self.x1, tf.random_normal(shape=tf.shape(self.x1), mean=0, stddev=0.2, dtype=tf.float32))
            x_input2 = tf.add(self.x2, tf.random_normal(shape=tf.shape(self.x2), mean=0, stddev=0.2, dtype=tf.float32))
            x_input3 = tf.add(self.x3, tf.random_normal(shape=tf.shape(self.x3), mean=0, stddev=0.2, dtype=tf.float32))

        latent1 = self.encoder(x_input1)
        latent2 = self.encoder(x_input2)
        latent3 = self.encoder(x_input3)

        self.z_conv1 = tf.reshape(latent1, [batch_size, -1])
        self.z_conv2 = tf.reshape(latent2, [batch_size, -1])
        self.z_conv3 = tf.reshape(latent3, [batch_size, -1])

        # self-expressive
        self.Coef_1 = tf.Variable(1.0e-8 * tf.ones([self.batch_size, self.batch_size], tf.float32), name='Coef_1')
        self.Coef_2 = tf.Variable(1.0e-8 * tf.ones([self.batch_size, self.batch_size], tf.float32), name='Coef_2')
        self.Coef_3 = tf.Variable(1.0e-8 * tf.ones([self.batch_size, self.batch_size], tf.float32), name='Coef_3')

        self.z1_c = tf.matmul(self.Coef_1, self.z_conv1)
        self.z2_c = tf.matmul(self.Coef_2, self.z_conv2)
        self.z3_c = tf.matmul(self.Coef_3, self.z_conv3)

        latent_de_ft1 = tf.reshape(self.z1_c, tf.shape(latent1))
        latent_de_ft2 = tf.reshape(self.z2_c, tf.shape(latent2))
        latent_de_ft3 = tf.reshape(self.z3_c, tf.shape(latent3))

        self.x_r1 = self.decoder(latent_de_ft1, x_input1)
        self.x_r2 = self.decoder(latent_de_ft2, x_input2)
        self.x_r3 = self.decoder(latent_de_ft3, x_input3)

        # Multi-view Fusion
        self.x4 = tf.concat((self.z_conv1, self.z_conv2, self.z_conv3), 1)
        # print(self.x4.shape)
        self.x4 = tf.reshape(self.x4, [batch_size, 1, self.x4.shape[1], 1])
        if denoise is False:
            x_input4 = self.x4
        else:
            x_input4 = tf.add(self.x4, tf.random_normal(shape=tf.shape(self.x4), mean=0, stddev=0.2, dtype=tf.float32))
        latent4 = self.encoder(x_input4)
        self.z_conv4 = tf.reshape(latent4, [batch_size, -1])
        self.Coef_4 = tf.Variable(1.0e-8 * tf.ones([self.batch_size, self.batch_size], tf.float32), name='Coef_4')
        self.z4_c = tf.matmul(self.Coef_4, self.z_conv4)
        latent_de_ft4 = tf.reshape(self.z4_c, tf.shape(latent4))
        self.x_r4 = self.decoder(latent_de_ft4, x_input4)

        self.saver = tf.train.Saver([v for v in tf.trainable_variables() if not (v.name.startswith("Coef"))])

        # self-expressive loss
        self.z_zc_loss1 = tf.reduce_sum(tf.pow(tf.subtract(self.z_conv1, self.z1_c), 2))
        self.z_zc_loss2 = tf.reduce_sum(tf.pow(tf.subtract(self.z_conv2, self.z2_c), 2))
        self.z_zc_loss3 = tf.reduce_sum(tf.pow(tf.subtract(self.z_conv3, self.z3_c), 2))
        self.z_zc_loss4 = tf.reduce_sum(tf.pow(tf.subtract(self.z_conv4, self.z4_c), 2))

        # recontruction loss
        self.recon_loss1 = 0.5 * tf.reduce_sum(tf.pow(tf.subtract(self.x_r1, x_input1), 2.0))
        self.recon_loss2 = 0.5 * tf.reduce_sum(tf.pow(tf.subtract(self.x_r2, x_input2), 2.0))
        self.recon_loss3 = 0.5 * tf.reduce_sum(tf.pow(tf.subtract(self.x_r3, x_input3), 2.0))
        self.recon_loss4 = 0.5 * tf.reduce_sum(tf.pow(tf.subtract(self.x_r4, x_input4), 2.0))

        # regular loss
        self.reg_loss1 = tf.reduce_sum(tf.pow(self.Coef_1, 2))
        self.reg_loss2 = tf.reduce_sum(tf.pow(self.Coef_2, 2))
        self.reg_loss3 = tf.reduce_sum(tf.pow(self.Coef_3, 2))
        self.reg_loss4 = tf.reduce_sum(tf.pow(self.Coef_4, 2))

        self.recon_loss = self.recon_loss1 + self.recon_loss2 + self.recon_loss3 + self.recon_loss4
        self.z_zc_loss = self.z_zc_loss1 + self.z_zc_loss2 + self.z_zc_loss3 + self.z_zc_loss4
        self.reg_loss = self.reg_loss1 + self.reg_loss2 + self.reg_loss3 + self.reg_loss4

        self.loss = self.recon_loss + reg_const1 * self.reg_loss + reg_const2 * self.z_zc_loss

        self.optimizer = tf.train.AdamOptimizer(learning_rate=self.learning_rate).minimize(self.loss)

        self.init = tf.global_variables_initializer()
        config = tf.ConfigProto()
        config.gpu_options.allow_growth = True
        self.sess = tf.InteractiveSession(config=config)
        self.sess.run(self.init)

    def _initialize_weights(self):
        all_weights = dict()

        all_weights['enc_w0'] = tf.Variable(para['W_encoder'])
        all_weights['enc_b0'] = tf.Variable(para['bias_e'])

        all_weights['dec_w0'] = tf.Variable(para['W_decoder'])
        all_weights['dec_b0'] = tf.Variable(para['bias_d'])

        return all_weights

    # Building the encoder
    def encoder(self, x):
        layer1 = tf.nn.relu(tf.nn.bias_add(tf.nn.conv2d(x, self.weights['enc_w0'], strides=[1, 1, 5, 1], padding='SAME'),
                            self.weights['enc_b0']))
        return layer1

    # Building the decoder
    def decoder(self, z, x):
        layer1 = tf.nn.relu(tf.add(tf.nn.conv2d_transpose(z, self.weights['dec_w0'], tf.stack(
                                               [tf.shape(x)[0], tf.shape(x)[1], tf.shape(x)[2], tf.shape(x)[3]]),
                                               strides=[1, 1, 5, 1], padding='SAME'), self.weights['dec_b0']))
        return layer1

    def finetune_fit(self, view1, view2, view3, lr):
        loss, Coef_1, Coef_2, Coef_3, Coef_4, _ = self.sess.run(
            (self.loss, self.Coef_1, self.Coef_2, self.Coef_3, self.Coef_4, self.optimizer),
            feed_dict={
                self.x1: view1,
                self.x2: view2,
                self.x3: view3,
                self.learning_rate: lr})
        return loss, Coef_1, Coef_2, Coef_3, Coef_4

    def initlization(self):
        # tf.reset_default_graph()
        self.sess.run(self.init)

    def get_latent(self, view1, view2, view3):
        latent_1, latent_2, latent_3 = self.sess.run(
            [self.z_conv1, self.z_conv2, self.z_conv3],
            feed_dict={
                self.x1: view1,
                self.x2: view2,
                self.x3: view3
            })
        return latent_1, latent_2, latent_3

    # def save_model(self):
    #     save_path = self.saver.save(self.sess, self.model_path)
    #     print("model_1 saved in file: %s" % save_path)
    #
    # def restore(self):
    #     self.saver.restore(self.sess, self.model_path)
    #     print("model_1 restored")


def thrC(C, ro):
    if ro < 1:
        N = C.shape[1]
        Cp = np.zeros((N, N))
        S = np.abs(np.sort(-np.abs(C), axis=0))
        Ind = np.argsort(-np.abs(C), axis=0)
        for i in range(N):
            cL1 = np.sum(S[:, i]).astype(float)
            stop = False
            csum = 0
            t = 0
            while (stop == False):
                csum = csum + S[t, i]
                if csum > ro * cL1:
                    stop = True
                    Cp[Ind[0:t + 1, i], i] = C[Ind[0:t + 1, i], i]
                t = t + 1
    else:
        Cp = C

    return Cp


def save_Coef(Coef, i, j):
    print(Coef)
    writer = pd.ExcelWriter(xlsx_dir + "model23/W%d_%d.xlsx" % (i, j))
    C = pd.DataFrame(Coef)
    C.to_excel(writer, 'page_1', float_format='%.8f')  # float_format 控制精度
    writer.save()
    print(xlsx_dir + 'model23/W%d_%d.xlsx' % (i, j))


def main():
    # 初始训练参数
    iter_ft = 0
    ft_times = 300
    # ft_times = 100  # tune
    display_step = 100
    alpha = 0.03  # tune
    learning_rate = 0.001  # tune

    reg1 = 1
    reg2 = 92
    kernel_size = [5]
    n_hidden = [15]

    for i in range(20):  # 总共训练多少次
        print("***************************************************")
        print("epoch:" + str(i))
        print("***************************************************")

        gene_expression0, Img0, Label0, n_input0 = load_data(data_name1)
        gene_expression1, Img1, Label1, n_input1 = load_data(data_name2)
        gene_expression2, Img2, Label2, n_input2 = load_data(data_name3)
        batch_size = gene_expression0.shape[0]

        CAE = ConvAE(n_input=[n_input0, n_input1, n_input2], n_hidden=n_hidden, reg_const1=reg1, reg_const2=reg2,
                     kernel_size=kernel_size, batch_size=batch_size)

        for _ in range(0, 1):
            CAE.initlization()
            # CAE.restore()

            for iter_ft in range(ft_times):
                iter_ft = iter_ft + 1
                recon_cost, Coef_1, Coef_2, Coef_3, Coef_4 = CAE.finetune_fit(Img0, Img1, Img2, lr=learning_rate)
                if iter_ft % display_step == 0:
                    print("epoch: %.1d" % iter_ft, "cost: %.8f" % (recon_cost / float(batch_size)))
                    C_1 = thrC(Coef_1, alpha)
                    C_2 = thrC(Coef_2, alpha)
                    C_3 = thrC(Coef_3, alpha)
                    C_4 = thrC(Coef_4, alpha)

        print("-" * 50)
        print("for the BREAST_Gene_Expression")
        save_Coef(C_1, i, 0)

        print("-" * 50)
        print("for the BREAST_Methy_Expression")
        save_Coef(C_2, i, 1)

        print("-" * 50)
        print("for the BREAST_Mirna_Expression")
        save_Coef(C_3, i, 2)

        print("-" * 50)
        print("for the Fusion_Expression")
        save_Coef(C_4, i, 3)

        file = open(xlsx_dir + 'model23/W%d.txt' % i, 'w')
        file.write("ft_times = " + str(ft_times) + " learning_rate = " + str(learning_rate) + " reg1 = " + str(
            reg1) + " reg2 = " + str(reg2) + " alpha = " + str(alpha))

        CAE.sess.close()

        # 下一次的训练参数
        iter_ft = 0
        # ft_times += 100
        ft_times += 300
        display_step = display_step
        alpha = alpha
        learning_rate = learning_rate

        reg1 = reg1
        # reg2 = reg2
        reg2 = reg2-1


if __name__ == "__main__":
    para = test_mnist(data_load(data_name1))  # 接着预训练再训练
    print("***************************************************")
    main()
