import tensorflow as tf
import os
import sys
import data_generation
import networks
import scipy.io as sio
import param
import util
import truncated_vgg
from keras.backend.tensorflow_backend import set_session
from keras.optimizers import Adam




def train(model_name, gpu_id):
    params = param.getGeneralParams()

    network_dir = params['project_dir'] + '/models/' + model_name

    if not os.path.isdir(network_dir):
        os.mkdir(network_dir)

    data_dir = '/afs/csail.mit.edu/u/b/balakg/pose/datasets/posewarp/train'
    train_feed = data_generation.create_feed(params, data_dir)

    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    set_session(tf.Session(config=config))

    vgg_model = truncated_vgg.vgg_norm()
    networks.make_trainable(vgg_model, False)
    response_weights = sio.loadmat('vgg_activation_distribution_train.mat')
    model = networks.network_posewarp(params)
    model.compile(optimizer=Adam(lr=1e-4), loss=[networks.vgg_loss(vgg_model, response_weights)])

    # model.summary()

    for step in xrange(0, 200000):
        x, y = next(train_feed)

        train_loss = model.train_on_batch(x, y)

        util.printProgress(step, 0, train_loss)

        if step > 0 and step % params['model_save_interval'] == 0:
            model.save(network_dir + '/' + str(step) + '.h5')


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "Need model name and gpu id as command line arguments."
    else:
        train(sys.argv[1], sys.argv[2])