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
/Fragments/chr9/RGLengths <HDF5 dataset "RGLengths": shape (1, 1717), type "<i4">
/Fragments/chr9/RGValues <HDF5 dataset "RGValues": shape (1, 1717), type "|S54">
/Fragments/chr9/Ranges <HDF5 dataset "Ranges": shape (2, 882526), type "<i4">
/GeneScoreMatrix/Info/CellNames <HDF5 dataset "CellNames": shape (1717,), type "|S54">
/GeneScoreMatrix/Info/Class <HDF5 dataset "Class": shape (1,), type "|S21">
/GeneScoreMatrix/Info/Date <HDF5 dataset "Date": shape (1,), type "|S11">
/GeneScoreMatrix/Info/FeatureDF <HDF5 dataset "FeatureDF": shape (27902,), type "|V44">
/GeneScoreMatrix/Info/Params <HDF5 dataset "Params": shape (2,), type "|V84">
/GeneScoreMatrix/Info/Units <HDF5 dataset "Units": shape (1,), type "|S11">
/GeneScoreMatrix/chr9/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/GeneScoreMatrix/chr9/i <HDF5 dataset "i": shape (1, 496633), type "<i4">
/GeneScoreMatrix/chr9/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr9/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/GeneScoreMatrix/chr9/rowMeansLog2 <HDF5 dataset "rowMeansLog2": shape (1, 1105), type "<f8">
/GeneScoreMatrix/chr9/rowSums <HDF5 dataset "rowSums": shape (1, 1105), type "<f8">
/GeneScoreMatrix/chr9/rowVarsLog2 <HDF5 dataset "rowVarsLog2": shape (1, 1105), type "<f8">
/GeneScoreMatrix/chr9/x <HDF5 dataset "x": shape (1, 496633), type "<f8">
/TileMatrix/chr9/colSums <HDF5 dataset "colSums": shape (1, 1717), type "<f8">
/TileMatrix/chr9/i <HDF5 dataset "i": shape (1, 1032893), type "<i4">
/TileMatrix/chr9/jLengths <HDF5 dataset "jLengths": shape (1, 1717), type "<i4">
/TileMatrix/chr9/jValues <HDF5 dataset "jValues": shape (1, 1717), type "<i4">
/TileMatrix/chr9/rowSums <HDF5 dataset "rowSums": shape (1, 266391), type "<f8">
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
```
