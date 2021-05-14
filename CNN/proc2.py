import warnings
warnings.filterwarnings('ignore')
import tensorflow as tf
FLAGS = tf.app.flags.FLAGS

from inception import nc_inception_eval
from inception.nc_imagenet_data import ImagenetData

f=open('s')
filenames = [line.rstrip() for line in f]
f.close()

def del_all_flags(FLAGS):
    flags_dict = FLAGS._flags()
    keys_list = [keys for keys in flags_dict]
    for keys in keys_list:
        FLAGS.__delattr__(keys)

del_all_flags(tf.app.flags.FLAGS)

path = '/home/wangqi9/tmp/t/'
tf.app.flags.DEFINE_string('checkpoint_dir', '/home/wangqi9/Downloads/checkpoints-20210416T004821Z-001/checkpoints/run1a_3D_classifier','')
tf.app.flags.DEFINE_string('eval_dir', path,'')
tf.app.flags.DEFINE_string('data_dir', path,'')
tf.app.flags.DEFINE_integer('batch_size', 100,'')
tf.app.flags.DEFINE_integer('image_size', 299,'')
tf.app.flags.DEFINE_integer('ClassNumber', 3,'')
tf.app.flags.DEFINE_string('mode', '0_softmax','')
tf.app.flags.DEFINE_string('TVmode', 'test','')
tf.app.flags.DEFINE_string('ImageSet_basename', 'test_','')
tf.app.flags.DEFINE_string('subset', 'valid', '')
tf.app.flags.DEFINE_integer('num_preprocess_threads', 4, '')
tf.app.flags.DEFINE_integer('num_readers', 4, '')
tf.app.flags.DEFINE_integer('num_examples', 1, '')
tf.app.flags.DEFINE_boolean('run_once', True, '')

dataset = ImagenetData(subset=FLAGS.subset)
dataset.data_files = lambda: filenames[15000:20000]
#dataset.data_files()

print('evaluating ...')
precision_at_1, current_score = nc_inception_eval.evaluate(dataset)
