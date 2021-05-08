## ArchR file format
```py
def h5py_dataset_iterator(g, prefix=''):
    for key in g.keys():
        item = g[key]
        path = '{}/{}'.format(prefix, key)
        if isinstance(item, h5py.Dataset): # test for dataset
            yield (path, item)
        elif isinstance(item, h5py.Group): # test for group (go down)
            yield from h5py_dataset_iterator(item, path)


f= h5py.File(s, 'r')
for (path, dset) in h5py_dataset_iterator(f):
	print(path, dset)
```

```
/ArchRVersion <HDF5 dataset "ArchRVersion": shape (1,), type "|S6">
/Class <HDF5 dataset "Class": shape (1,), type "|S6">
/Fragments/chr1/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chr1/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chr1/Ranges <HDF5 dataset "Ranges": shape (2, 1780973), type "<i4">
/Fragments/chr10/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chr10/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chr10/Ranges <HDF5 dataset "Ranges": shape (2, 913212), type "<i4">
/Fragments/chr11/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chr11/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chr11/Ranges <HDF5 dataset "Ranges": shape (2, 1008393), type "<i4">
/Fragments/chr12/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chr12/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chr12/Ranges <HDF5 dataset "Ranges": shape (2, 745640), type "<i4">
/Fragments/chr13/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chr13/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chr13/Ranges <HDF5 dataset "Ranges": shape (2, 771124), type "<i4">
/Fragments/chr14/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chr14/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chr14/Ranges <HDF5 dataset "Ranges": shape (2, 1062279), type "<i4">
/Fragments/chr15/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chr15/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chr15/Ranges <HDF5 dataset "Ranges": shape (2, 813376), type "<i4">
/Fragments/chr16/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chr16/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chr16/Ranges <HDF5 dataset "Ranges": shape (2, 890234), type "<i4">
/Fragments/chr17/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chr17/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chr17/Ranges <HDF5 dataset "Ranges": shape (2, 491537), type "<i4">
/Fragments/chr18/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chr18/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chr18/Ranges <HDF5 dataset "Ranges": shape (2, 418786), type "<i4">
/Fragments/chr19/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chr19/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chr19/Ranges <HDF5 dataset "Ranges": shape (2, 788615), type "<i4">
/Fragments/chr2/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chr2/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chr2/Ranges <HDF5 dataset "Ranges": shape (2, 1184799), type "<i4">
/Fragments/chr20/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chr20/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chr20/Ranges <HDF5 dataset "Ranges": shape (2, 730891), type "<i4">
/Fragments/chr3/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chr3/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chr3/Ranges <HDF5 dataset "Ranges": shape (2, 1174385), type "<i4">
/Fragments/chr4/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chr4/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chr4/Ranges <HDF5 dataset "Ranges": shape (2, 1008266), type "<i4">
/Fragments/chr5/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chr5/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chr5/Ranges <HDF5 dataset "Ranges": shape (2, 872400), type "<i4">
/Fragments/chr6/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chr6/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chr6/Ranges <HDF5 dataset "Ranges": shape (2, 1099990), type "<i4">
/Fragments/chr7/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chr7/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chr7/Ranges <HDF5 dataset "Ranges": shape (2, 1236054), type "<i4">
/Fragments/chr8/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chr8/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chr8/Ranges <HDF5 dataset "Ranges": shape (2, 919517), type "<i4">
/Fragments/chr9/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chr9/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chr9/Ranges <HDF5 dataset "Ranges": shape (2, 882526), type "<i4">
/Fragments/chrX/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chrX/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chrX/Ranges <HDF5 dataset "Ranges": shape (2, 483374), type "<i4">
/GeneScoreMatrix/Info/CellNames <HDF5 dataset "CellNames": shape (1717,), type "|S54">
/GeneScoreMatrix/Info/Class <HDF5 dataset "Class": shape (1,), type "|S21">
/GeneScoreMatrix/Info/Date <HDF5 dataset "Date": shape (1,), type "|S11">
/GeneScoreMatrix/Info/FeatureDF <HDF5 dataset "FeatureDF": shape (27902,), type "|V44">
/GeneScoreMatrix/Info/Params <HDF5 dataset "Params": shape (2,), type "|V84">
/GeneScoreMatrix/Info/Units <HDF5 dataset "Units": shape (1,), type "|S11">
/GeneScoreMatrix/chr1/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chr1/i <HDF5 dataset "i": shape (1, 1106093), type "<i4">
/GeneScoreMatrix/chr1/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr1/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr1/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 2737), type "<f8">
/GeneScoreMatrix/chr1/rowSums <HDF5 dataset "rowSums": shape (1, 2737), type "<f8">
/GeneScoreMatrix/chr1/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 2737), type "<f8">
/GeneScoreMatrix/chr1/x <HDF5 dataset "x": shape (1, 1106093), type "<f8">
/GeneScoreMatrix/chr10/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chr10/i <HDF5 dataset "i": shape (1, 577126), type "<i4">
/GeneScoreMatrix/chr10/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr10/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr10/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 1379), type "<f8">
/GeneScoreMatrix/chr10/rowSums <HDF5 dataset "rowSums": shape (1, 1379), type "<f8">
/GeneScoreMatrix/chr10/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 1379), type "<f8">
/GeneScoreMatrix/chr10/x <HDF5 dataset "x": shape (1, 577126), type "<f8">
/GeneScoreMatrix/chr11/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chr11/i <HDF5 dataset "i": shape (1, 615821), type "<i4">
/GeneScoreMatrix/chr11/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr11/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr11/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 1498), type "<f8">
/GeneScoreMatrix/chr11/rowSums <HDF5 dataset "rowSums": shape (1, 1498), type "<f8">
/GeneScoreMatrix/chr11/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 1498), type "<f8">
/GeneScoreMatrix/chr11/x <HDF5 dataset "x": shape (1, 615821), type "<f8">
/GeneScoreMatrix/chr12/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chr12/i <HDF5 dataset "i": shape (1, 434759), type "<i4">
/GeneScoreMatrix/chr12/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr12/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr12/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 992), type "<f8">
/GeneScoreMatrix/chr12/rowSums <HDF5 dataset "rowSums": shape (1, 992), type "<f8">
/GeneScoreMatrix/chr12/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 992), type "<f8">
/GeneScoreMatrix/chr12/x <HDF5 dataset "x": shape (1, 434759), type "<f8">
/GeneScoreMatrix/chr13/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chr13/i <HDF5 dataset "i": shape (1, 427486), type "<i4">
/GeneScoreMatrix/chr13/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr13/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr13/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 987), type "<f8">
/GeneScoreMatrix/chr13/rowSums <HDF5 dataset "rowSums": shape (1, 987), type "<f8">
/GeneScoreMatrix/chr13/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 987), type "<f8">
/GeneScoreMatrix/chr13/x <HDF5 dataset "x": shape (1, 427486), type "<f8">
/GeneScoreMatrix/chr14/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chr14/i <HDF5 dataset "i": shape (1, 653705), type "<i4">
/GeneScoreMatrix/chr14/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr14/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr14/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 1580), type "<f8">
/GeneScoreMatrix/chr14/rowSums <HDF5 dataset "rowSums": shape (1, 1580), type "<f8">
/GeneScoreMatrix/chr14/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 1580), type "<f8">
/GeneScoreMatrix/chr14/x <HDF5 dataset "x": shape (1, 653705), type "<f8">
/GeneScoreMatrix/chr15/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chr15/i <HDF5 dataset "i": shape (1, 469284), type "<i4">
/GeneScoreMatrix/chr15/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr15/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr15/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 1049), type "<f8">
/GeneScoreMatrix/chr15/rowSums <HDF5 dataset "rowSums": shape (1, 1049), type "<f8">
/GeneScoreMatrix/chr15/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 1049), type "<f8">
/GeneScoreMatrix/chr15/x <HDF5 dataset "x": shape (1, 469284), type "<f8">
/GeneScoreMatrix/chr16/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chr16/i <HDF5 dataset "i": shape (1, 605461), type "<i4">
/GeneScoreMatrix/chr16/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr16/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr16/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 1486), type "<f8">
/GeneScoreMatrix/chr16/rowSums <HDF5 dataset "rowSums": shape (1, 1486), type "<f8">
/GeneScoreMatrix/chr16/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 1486), type "<f8">
/GeneScoreMatrix/chr16/x <HDF5 dataset "x": shape (1, 605461), type "<f8">
/GeneScoreMatrix/chr17/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chr17/i <HDF5 dataset "i": shape (1, 247500), type "<i4">
/GeneScoreMatrix/chr17/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr17/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr17/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 546), type "<f8">
/GeneScoreMatrix/chr17/rowSums <HDF5 dataset "rowSums": shape (1, 546), type "<f8">
/GeneScoreMatrix/chr17/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 546), type "<f8">
/GeneScoreMatrix/chr17/x <HDF5 dataset "x": shape (1, 247500), type "<f8">
/GeneScoreMatrix/chr18/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chr18/i <HDF5 dataset "i": shape (1, 211593), type "<i4">
/GeneScoreMatrix/chr18/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr18/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr18/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 470), type "<f8">
/GeneScoreMatrix/chr18/rowSums <HDF5 dataset "rowSums": shape (1, 470), type "<f8">
/GeneScoreMatrix/chr18/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 470), type "<f8">
/GeneScoreMatrix/chr18/x <HDF5 dataset "x": shape (1, 211593), type "<f8">
/GeneScoreMatrix/chr19/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chr19/i <HDF5 dataset "i": shape (1, 592160), type "<i4">
/GeneScoreMatrix/chr19/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr19/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr19/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 1695), type "<f8">
/GeneScoreMatrix/chr19/rowSums <HDF5 dataset "rowSums": shape (1, 1695), type "<f8">
/GeneScoreMatrix/chr19/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 1695), type "<f8">
/GeneScoreMatrix/chr19/x <HDF5 dataset "x": shape (1, 592160), type "<f8">
/GeneScoreMatrix/chr2/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chr2/i <HDF5 dataset "i": shape (1, 656291), type "<i4">
/GeneScoreMatrix/chr2/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr2/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr2/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 1545), type "<f8">
/GeneScoreMatrix/chr2/rowSums <HDF5 dataset "rowSums": shape (1, 1545), type "<f8">
/GeneScoreMatrix/chr2/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 1545), type "<f8">
/GeneScoreMatrix/chr2/x <HDF5 dataset "x": shape (1, 656291), type "<f8">
/GeneScoreMatrix/chr20/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chr20/i <HDF5 dataset "i": shape (1, 438647), type "<i4">
/GeneScoreMatrix/chr20/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr20/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr20/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 1021), type "<f8">
/GeneScoreMatrix/chr20/rowSums <HDF5 dataset "rowSums": shape (1, 1021), type "<f8">
/GeneScoreMatrix/chr20/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 1021), type "<f8">
/GeneScoreMatrix/chr20/x <HDF5 dataset "x": shape (1, 438647), type "<f8">
/GeneScoreMatrix/chr3/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chr3/i <HDF5 dataset "i": shape (1, 674003), type "<i4">
/GeneScoreMatrix/chr3/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr3/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr3/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 1585), type "<f8">
/GeneScoreMatrix/chr3/rowSums <HDF5 dataset "rowSums": shape (1, 1585), type "<f8">
/GeneScoreMatrix/chr3/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 1585), type "<f8">
/GeneScoreMatrix/chr3/x <HDF5 dataset "x": shape (1, 674003), type "<f8">
/GeneScoreMatrix/chr4/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chr4/i <HDF5 dataset "i": shape (1, 600282), type "<i4">
/GeneScoreMatrix/chr4/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr4/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr4/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 1468), type "<f8">
/GeneScoreMatrix/chr4/rowSums <HDF5 dataset "rowSums": shape (1, 1468), type "<f8">
/GeneScoreMatrix/chr4/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 1468), type "<f8">
/GeneScoreMatrix/chr4/x <HDF5 dataset "x": shape (1, 600282), type "<f8">
/GeneScoreMatrix/chr5/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chr5/i <HDF5 dataset "i": shape (1, 467075), type "<i4">
/GeneScoreMatrix/chr5/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr5/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr5/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 1182), type "<f8">
/GeneScoreMatrix/chr5/rowSums <HDF5 dataset "rowSums": shape (1, 1182), type "<f8">
/GeneScoreMatrix/chr5/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 1182), type "<f8">
/GeneScoreMatrix/chr5/x <HDF5 dataset "x": shape (1, 467075), type "<f8">
/GeneScoreMatrix/chr6/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chr6/i <HDF5 dataset "i": shape (1, 583408), type "<i4">
/GeneScoreMatrix/chr6/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr6/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr6/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 1318), type "<f8">
/GeneScoreMatrix/chr6/rowSums <HDF5 dataset "rowSums": shape (1, 1318), type "<f8">
/GeneScoreMatrix/chr6/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 1318), type "<f8">
/GeneScoreMatrix/chr6/x <HDF5 dataset "x": shape (1, 583408), type "<f8">
/GeneScoreMatrix/chr7/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chr7/i <HDF5 dataset "i": shape (1, 758960), type "<i4">
/GeneScoreMatrix/chr7/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr7/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr7/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 1953), type "<f8">
/GeneScoreMatrix/chr7/rowSums <HDF5 dataset "rowSums": shape (1, 1953), type "<f8">
/GeneScoreMatrix/chr7/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 1953), type "<f8">
/GeneScoreMatrix/chr7/x <HDF5 dataset "x": shape (1, 758960), type "<f8">
/GeneScoreMatrix/chr8/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chr8/i <HDF5 dataset "i": shape (1, 474428), type "<i4">
/GeneScoreMatrix/chr8/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr8/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr8/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 1036), type "<f8">
/GeneScoreMatrix/chr8/rowSums <HDF5 dataset "rowSums": shape (1, 1036), type "<f8">
/GeneScoreMatrix/chr8/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 1036), type "<f8">
/GeneScoreMatrix/chr8/x <HDF5 dataset "x": shape (1, 474428), type "<f8">
/GeneScoreMatrix/chr9/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chr9/i <HDF5 dataset "i": shape (1, 496633), type "<i4">
/GeneScoreMatrix/chr9/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr9/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr9/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 1105), type "<f8">
/GeneScoreMatrix/chr9/rowSums <HDF5 dataset "rowSums": shape (1, 1105), type "<f8">
/GeneScoreMatrix/chr9/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 1105), type "<f8">
/GeneScoreMatrix/chr9/x <HDF5 dataset "x": shape (1, 496633), type "<f8">
/GeneScoreMatrix/chrX/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chrX/i <HDF5 dataset "i": shape (1, 318170), type "<i4">
/GeneScoreMatrix/chrX/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chrX/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chrX/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 1270), type "<f8">
/GeneScoreMatrix/chrX/rowSums <HDF5 dataset "rowSums": shape (1, 1270), type "<f8">
/GeneScoreMatrix/chrX/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 1270), type "<f8">
/GeneScoreMatrix/chrX/x <HDF5 dataset "x": shape (1, 318170), type "<f8">
/Metadata/BlacklistRatio <HDF5 dataset "BlacklistRatio": shape (1717,), type "<f8">
/Metadata/CellNames <HDF5 dataset "CellNames": shape (1717,), type "|S54">
/Metadata/Completed <HDF5 dataset "Completed": shape (1,), type "|S9">
/Metadata/Date <HDF5 dataset "Date": shape (1,), type "|S11">
/Metadata/DoubletEnrichment <HDF5 dataset "DoubletEnrichment": shape (1717,), type "<f8">
/Metadata/DoubletScore <HDF5 dataset "DoubletScore": shape (1717,), type "<f8">
/Metadata/NucleosomeRatio <HDF5 dataset "NucleosomeRatio": shape (1717,), type "<f8">
/Metadata/PassQC <HDF5 dataset "PassQC": shape (1717,), type "<f8">
/Metadata/PromoterRatio <HDF5 dataset "PromoterRatio": shape (1717,), type "<f8">
/Metadata/ReadsInBlacklist <HDF5 dataset "ReadsInBlacklist": shape (1717,), type "<f8">
/Metadata/ReadsInPromoter <HDF5 dataset "ReadsInPromoter": shape (1717,), type "<f8">
/Metadata/ReadsInTSS <HDF5 dataset "ReadsInTSS": shape (1717,), type "<f8">
/Metadata/Sample <HDF5 dataset "Sample": shape (1,), type "|S43">
/Metadata/TSSEnrichment <HDF5 dataset "TSSEnrichment": shape (1717,), type "<f8">
/Metadata/nDiFrags <HDF5 dataset "nDiFrags": shape (1717,), type "<f8">
/Metadata/nFrags <HDF5 dataset "nFrags": shape (1717,), type "<f8">
/Metadata/nMonoFrags <HDF5 dataset "nMonoFrags": shape (1717,), type "<f8">
/Metadata/nMultiFrags <HDF5 dataset "nMultiFrags": shape (1717,), type "<f8">
/TileMatrix/Info/CellNames <HDF5 dataset "CellNames": shape (1717,), type "|S54">
/TileMatrix/Info/Class <HDF5 dataset "Class": shape (1,), type "|S21">
/TileMatrix/Info/Date <HDF5 dataset "Date": shape (1,), type "|S11">
/TileMatrix/Info/FeatureDF <HDF5 dataset "FeatureDF": shape (5743663,), type "|V18">
/TileMatrix/Info/Params <HDF5 dataset "Params": shape (21,), type "|V18">
/TileMatrix/Info/Units <HDF5 dataset "Units": shape (1,), type "|S16">
/TileMatrix/chr1/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chr1/i <HDF5 dataset "i": shape (1, 2066363), type "<i4">
/TileMatrix/chr1/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chr1/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chr1/rowSums <HDF5 dataset "rowSums": shape (1, 455113), type "<f8">
/TileMatrix/chr10/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chr10/i <HDF5 dataset "i": shape (1, 1051323), type "<i4">
/TileMatrix/chr10/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chr10/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chr10/rowSums <HDF5 dataset "rowSums": shape (1, 193020), type "<f8">
/TileMatrix/chr11/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chr11/i <HDF5 dataset "i": shape (1, 1179336), type "<i4">
/TileMatrix/chr11/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chr11/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chr11/rowSums <HDF5 dataset "rowSums": shape (1, 275516), type "<f8">
/TileMatrix/chr12/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chr12/i <HDF5 dataset "i": shape (1, 881486), type "<i4">
/TileMatrix/chr12/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chr12/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chr12/rowSums <HDF5 dataset "rowSums": shape (1, 265174), type "<f8">
/TileMatrix/chr13/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chr13/i <HDF5 dataset "i": shape (1, 899480), type "<i4">
/TileMatrix/chr13/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chr13/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chr13/rowSums <HDF5 dataset "rowSums": shape (1, 222387), type "<f8">
/TileMatrix/chr14/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chr14/i <HDF5 dataset "i": shape (1, 1233034), type "<i4">
/TileMatrix/chr14/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chr14/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chr14/rowSums <HDF5 dataset "rowSums": shape (1, 261467), type "<f8">
/TileMatrix/chr15/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chr15/i <HDF5 dataset "i": shape (1, 947264), type "<i4">
/TileMatrix/chr15/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chr15/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chr15/rowSums <HDF5 dataset "rowSums": shape (1, 225226), type "<f8">
/TileMatrix/chr16/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chr16/i <HDF5 dataset "i": shape (1, 1014408), type "<i4">
/TileMatrix/chr16/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chr16/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chr16/rowSums <HDF5 dataset "rowSums": shape (1, 161996), type "<f8">
/TileMatrix/chr17/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chr17/i <HDF5 dataset "i": shape (1, 585453), type "<i4">
/TileMatrix/chr17/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chr17/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chr17/rowSums <HDF5 dataset "rowSums": shape (1, 193730), type "<f8">
/TileMatrix/chr18/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chr18/i <HDF5 dataset "i": shape (1, 497263), type "<i4">
/TileMatrix/chr18/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chr18/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chr18/rowSums <HDF5 dataset "rowSums": shape (1, 151424), type "<f8">
/TileMatrix/chr19/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chr19/i <HDF5 dataset "i": shape (1, 876792), type "<i4">
/TileMatrix/chr19/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chr19/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chr19/rowSums <HDF5 dataset "rowSums": shape (1, 118497), type "<f8">
/TileMatrix/chr2/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chr2/i <HDF5 dataset "i": shape (1, 1388942), type "<i4">
/TileMatrix/chr2/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chr2/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chr2/rowSums <HDF5 dataset "rowSums": shape (1, 384921), type "<f8">
/TileMatrix/chr20/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chr20/i <HDF5 dataset "i": shape (1, 819459), type "<i4">
/TileMatrix/chr20/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chr20/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chr20/rowSums <HDF5 dataset "rowSums": shape (1, 157083), type "<f8">
/TileMatrix/chr3/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chr3/i <HDF5 dataset "i": shape (1, 1377506), type "<i4">
/TileMatrix/chr3/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chr3/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chr3/rowSums <HDF5 dataset "rowSums": shape (1, 384589), type "<f8">
/TileMatrix/chr4/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chr4/i <HDF5 dataset "i": shape (1, 1194397), type "<i4">
/TileMatrix/chr4/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chr4/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chr4/rowSums <HDF5 dataset "rowSums": shape (1, 341911), type "<f8">
/TileMatrix/chr5/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chr5/i <HDF5 dataset "i": shape (1, 1037171), type "<i4">
/TileMatrix/chr5/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chr5/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chr5/rowSums <HDF5 dataset "rowSums": shape (1, 378909), type "<f8">
/TileMatrix/chr6/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chr6/i <HDF5 dataset "i": shape (1, 1301253), type "<i4">
/TileMatrix/chr6/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chr6/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chr6/rowSums <HDF5 dataset "rowSums": shape (1, 363170), type "<f8">
/TileMatrix/chr7/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chr7/i <HDF5 dataset "i": shape (1, 1440206), type "<i4">
/TileMatrix/chr7/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chr7/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chr7/rowSums <HDF5 dataset "rowSums": shape (1, 343765), type "<f8">
/TileMatrix/chr8/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chr8/i <HDF5 dataset "i": shape (1, 1062971), type "<i4">
/TileMatrix/chr8/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chr8/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chr8/rowSums <HDF5 dataset "rowSums": shape (1, 293702), type "<f8">
/TileMatrix/chr9/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chr9/i <HDF5 dataset "i": shape (1, 1032893), type "<i4">
/TileMatrix/chr9/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chr9/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chr9/rowSums <HDF5 dataset "rowSums": shape (1, 266391), type "<f8">
/TileMatrix/chrX/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chrX/i <HDF5 dataset "i": shape (1, 572102), type "<i4">
/TileMatrix/chrX/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chrX/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chrX/rowSums <HDF5 dataset "rowSums": shape (1, 305672), type "<f8">
```
