import tensorflow as tf
import os,sys

from build_TF_test import ImageCoder, _process_image, _convert_to_example

#filenames = [f for f in os.listdir('./') if f.endswith('.jpeg')]
f=open('2')
filenames = [line.rstrip() for line in f]
f.close()
coder = ImageCoder()

for filename in filenames:
    image_buffer, height, width = _process_image( filename, coder)
    example = _convert_to_example(filename, image_buffer, 1, 'r', height, width)
    writer = tf.python_io.TFRecordWriter(os.path.splitext(filename)[0]+'.TFRecord')
    writer.write(example.SerializeToString())
    writer.close()
