import tensorflow as tf
FLAGS = tf.app.flags.FLAGS

from inception import nc_inception_eval
from inception.nc_imagenet_data import ImagenetData
import numpy as np
import os,sys

def del_all_flags(FLAGS):
    flags_dict = FLAGS._flags()
    keys_list = [keys for keys in flags_dict]
    for keys in keys_list:
        FLAGS.__delattr__(keys)
        
from build_TF_test import ImageCoder, _process_image, _convert_to_example

path=sys.argv[1]
filenames = [f for f in os.listdir(path) if f.endswith('.jpeg')]
coder = ImageCoder()

for filename in filenames:
    image_buffer, height, width = _process_image(path + filename, coder)
    example = _convert_to_example(filename, image_buffer, 1, 'r', height, width)
    writer = tf.python_io.TFRecordWriter(path+'test_'+os.path.splitext(filename)[0]+'.TFRecord')
    writer.write(example.SerializeToString())
    writer.close()

sys.exit(0)
del_all_flags(FLAGS)

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
tf.app.flags.DEFINE_boolean('nbr_of_classes', False, '')

dataset = ImagenetData(subset=FLAGS.subset)
dataset.data_files = lambda: tf.gfile.Glob(os.path.join(path, '*.TFRecord'))
#dataset.data_files()

precision_at_1, current_score = nc_inception_eval.evaluate(dataset)

