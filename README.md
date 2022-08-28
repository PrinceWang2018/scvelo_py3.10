# scVelo在python3.10环境下的安装及使用

最近在处理单细胞测序的数据时需要分析RNA速率的变化，但在处理时遇到了一些安装问题跟大家分享。

走常规单细胞测序处理流程，跑完CellRanger，得到输出文件夹及比对后的bam文件，按照常规的方法合并、降维、注释并取亚群。

分析RNA velocity主要用的是velocyto.py的python程序，这个程序很好安装，参考官网流程即可。

> http://velocyto.org/velocyto.py/install/index.html

``` shell
#在创建的conda环境中（python 3.10也可以）
conda install numpy scipy cython numba matplotlib scikit-learn h5py click
pip install pysam
pip install velocyto
```

## 文件准备

1. cellranger运行后的输出文件
2. 标准的基因组gtf注释文件，可以在cellranger的reference文件中找到：refdata-gex-GRCh38-2020-A/genes/genes.gtf
3. 对repetetive elements注释的gtf文件 ，可以从UCSC官网下载，选择对应物种版本，track选择RepeatMasker，table选择rmsk，输出gtf文件：hg38_repeat_rmsk.gtf

> https://genome.ucsc.edu/cgi-bin/hgTables

## 数据预处理

需将possorted_genome_bam.bam按照cell barcode （-t CB）重新排序，生成cellsorted_possorted_genome_bam.bam，并放在相同目录下

```shell
samtools sort -@ 8  -t CB -O BAM -o cellsorted_possorted_genome_bam.bam possorted_genome_bam.bam
```

## 运行velocyto

```shell
#在安装velocyto的环境下
cellranger_outDir=/cellranger_outDir/   #不要具体到out文件夹
cellranger_gtf=/cellranger_ref_dir/refdata-gex-GRCh38-2020-A/genes/genes.gtf
rmsk_gtf=/rmsk_gtf_dir/hg38_repeat_rmsk.gtf
ls -lh $rmsk_gtf  $cellranger_outDir $cellranger_gtf
velocyto run10x -m $rmsk_gtf  $cellranger_outDir $cellranger_gtf
#这一步比较耗时
```

正常情况下会在漫长的比对后输出.loom文件，大小约100-200M。

## 使用loom文件进行下一步的处理-velocyto.R

下一步的分析最开始了解到的还是velocyto.R的R包，据说这个包很难装，最后尝试使用conda安装成功了。但是看到出图及分析也很有限，还是没有用他继续分析。

```
conda install -c bioconda r-velocyto.r
```

在实际使用过程中还需要同时使用Seurat及SeuratWrappers等众多R包，从头安装势必会遇到一点点挑战。

## scVelo的安装

而另一个分析RNA速率的python程序scVelo分析更全面，出图也很漂亮，教程也都很简单，但安装时候还是遇到了麻烦。

> https://scvelo.readthedocs.io/installation/

官网给出的使用教程中给出了pip的安装方式：

```sh
pip install -U scvelo
# scVelo requires Python 3.6 or later.
```

但实际上在python 3.10环境中即使安装成功了，import时候也会有非常多的报错，基本都是版本问题，试着解决了几个还是放弃了。后来尝试conda安装，果然还是conda大法好，不知道官网为什么没有推荐，这里使用的也是默认的3.10版本的python。

```shell
conda install scvelo   #默认0.2.24版本的scelo
#同时安装了anndata、pyscan、loompy等python包的依赖
python
```

进入python

``` python
import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt

#导入anndata时候会报错
>>> sample=anndata.read_loom("/loomdir/Sample.loom")
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/home/wzx/miniconda3/envs/velocyto/lib/python3.10/site-packages/anndata/compat/__init__.py", line 305, in inner_f
    return f(*args, **kwargs)
  File "/home/wzx/miniconda3/envs/velocyto/lib/python3.10/site-packages/anndata/_io/read.py", line 251, in read_loom
    with connect(filename, "r", **kwargs) as lc:
  File "/home/wzx/miniconda3/envs/velocyto/lib/python3.10/site-packages/loompy/loompy.py", line 1152, in connect
    return LoomConnection(filename, mode, validate=validate, spec_version=spec_version)
  File "/home/wzx/miniconda3/envs/velocyto/lib/python3.10/site-packages/loompy/loompy.py", line 83, in __init__
    raise ValueError("\n".join(lv.errors) + f"\n{filename} does not appead to be a valid Loom file according to Loom spec version '{spec_version}'")
ValueError: Row attribute 'Accession' dtype object is not allowed
Row attribute 'Chromosome' dtype object is not allowed
Row attribute 'Gene' dtype object is not allowed
Row attribute 'Strand' dtype object is not allowed
Column attribute 'CellID' dtype object is not allowed
For help, see http://linnarssonlab.org/loompy/format/
Object.loom does not appead to be a valid Loom file according to Loom spec version '2.0.1'
```

最开始怀疑是自己做的loom文件有问题，但反复检查了velocity的log文件，除了几个被允许的warning基本没有error。loom文件也不能直接打开编辑，cat会出现乱码。从GEO数据库下载了几份别人的loom文件发现同样报错，说明是scVelo安装的问题。如报错所说，可能是loompy.py软件内部的错误，velocity生成的loom文件不符合最新的loompy导入的标准，大概查了以下，导入前是要进行文件检查：

> Loom files (`.loom`) are created in the [HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format) file format, which supports an internal collection of numerical multidimensional datasets. HDF5 is supported by many computer languages, including [Python](http://h5py.org/), [R](http://bioconductor.org/packages/release/bioc/html/rhdf5.html), [MATLAB](http://se.mathworks.com/help/matlab/low-level-functions.html), [Mathematica](https://reference.wolfram.com/language/ref/format/HDF5.html), [C](https://www.hdfgroup.org/HDF5/doc/index.html), [C++](https://www.hdfgroup.org/HDF5/doc/cpplus_RM/), [Java](https://www.hdfgroup.org/products/java/), and [Ruby](https://rubygems.org/gems/hdf5/versions/0.3.5).
>
> Loom v3.0.0 introduces two major backwards-incompatible changes (global attributes and variable-length strings; see above).
>
> A compliant Loom reader MUST check the LOOM_SPEC_VERSION and treat files consistently with their spec. For example, when writing a global attribute, the writer MUST write only to the `/attrs` group if `LOOM_SPEC_VERSION` is `3.0.0` or higher. The writer MUST write the HDF5 attributes on the root `/` group if `LOOM_SPEC_VERSION` is lower than `3.0.0` or if it does not exist. This is to preserve a consistent format for legacy files.

既然是新标准，那么肯定有对应的旧版本软件能满足正常运行的要求，所以尝试着对loompy降级，conda默认安装的loompy是2.0.16，逐渐安装更低级别的包，随机试了2.0.15，2.0.14，2.0.6几个版本，在2.0.6时报错不兼容：

``` shell
ERROR: pip's dependency resolver does not currently take into account all the packages that are installed. This behaviour is the source of the following dependency conflicts.
scvelo 0.2.4 requires loompy>=2.0.12, but you have loompy 2.0.6 which is incompatible.
```

提示scvelo要求loompy最低版本2.0.12，安装后再次导入loom文件的时候则成功！！！

``` shell
pip install loompy==2.0.12
```

## 从Seurat Object提取坐标及cluster信息-R

运行流程部分参考@Yu_Lin老师的教程，感谢老师的分享：

> https://www.jianshu.com/p/c33341b65cad

```r
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(Seurat)
End.combined<-readRDS("./RDS.Rdata")
Eutopic<-End.combined

#修改Cells ID与loom文件一致
df<-data.frame(Cells=Cells(End.combined),sample=paste(End.combined@meta.data[["sample"]])) 
df$id<-sapply(df$Cells,function(x)paste(unlist(strsplit(x, "-"))[1],"x",sep = ""))
df$Cells<-paste(df$sample,df$id,sep = ":")
write.csv(df$Cells, file = "cellID_obs.csv", row.names = FALSE) 
              
## 提取坐标及相关信息
cell_embeddings<-Embeddings(Eutopic, reduction = "umap")
rownames(cell_embeddings)<-df$Cells 
write.csv(cell_embeddings, file = "cell_embeddings.csv")
clusters_obs<-Eutopic$seurat_subclusters 
names(clusters_obs)<-df$Cells 
write.csv(clusters_obs, file = "clusters_obs.csv")
```

## 运行scVelo

``` python
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt

# Import loom files
sample = anndata.read_loom("/loomdir/sample.loom")

# Import Seurat data
sample_obs = pd.read_csv("/cellID_dir/cellID_obs.csv")
umap = pd.read_csv("/cell_embeddings_dir/cell_embeddings.csv")
cell_clusters = pd.read_csv("/clusters_dir/clusters_obs.csv")

# Filter anndata according to Seurat data
sample = sample[np.isin(sample.obs.index,sample_obs["x"])]

#Unique variable names
sample.var_names_make_unique()

#Add Seurat metadata to anndata
adata_index = pd.DataFrame(adata.obs.index)
adata_index = adata_index.rename(columns = {0:'Cell ID'})
adata_index = adata_index.rename(columns = {"CellID":'Cell ID'})
rep=lambda x : x.split("-")[0]
adata_index["Cell ID"]=adata_index["Cell ID"].apply(rep)

#Add UMAP data
umap = umap.rename(columns = {'Unnamed: 0':'Cell ID'}) 
umap = umap[np.isin(umap["Cell ID"],adata_index["Cell ID"])]   
umap=umap.drop_duplicates(subset=["Cell ID"])  
umap_ordered = adata_index.merge(umap, on = "Cell ID")  
umap_ordered = umap_ordered.iloc[:,1:]  
adata.obsm['X_umap'] = umap_ordered.values   

#Add cluster data
cell_clusters = cell_clusters.rename(columns = {'Unnamed: 0':'Cell ID'}) 
cell_clusters = cell_clusters[np.isin(cell_clusters["Cell ID"],adata_index["Cell ID"])]
cell_clusters=cell_clusters.drop_duplicates(subset=["Cell ID"]) 
cell_clusters_ordered = adata_index.merge(cell_clusters, on = "Cell ID")
cell_clusters_ordered = cell_clusters_ordered.iloc[:,1:]
adata.obs['clusters']=cell_clusters_ordered.values
```

后续分析参考官网说明，放一部分代码供大家参考：

> https://scvelo.readthedocs.io/VelocityBasics/

``` python
# Basic processing
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
scv.tl.velocity(adata, mode = "stochastic")
scv.tl.velocity_graph(adata)

# spliced/unspliced的比例
scv.pl.proportions(adata)

# Visualization
scv.pl.velocity_embedding(adata, basis = 'umap')
scv.pl.velocity_embedding_stream(adata, basis = 'umap')
scv.pl.velocity(adata, ['geneA','geneB','geneC','geneD'], ncols=2)
scv.pl.scatter(adata, 'Cpe', color=['clusters', 'velocity'],
               add_outline='Ngn3 high EP, Pre-endocrine, Beta')

# Identify rank_velocity_genes
scv.tl.rank_velocity_genes(adata, groupby='clusters', min_corr=.3)
df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()
df.to_csv('/output/Rank_velocity_genes.csv', index = False)

kwargs = dict(frameon=False, size=10, linewidth=1.5,
              add_outline='Ngn3 high EP, Pre-endocrine, Beta')
scv.pl.scatter(adata, df['Ngn3 high EP'][:5], ylabel='Ngn3 high EP', **kwargs)
scv.pl.scatter(adata, df['Pre-endocrine'][:5], ylabel='Pre-endocrine', **kwargs)

#Velocities in cycling progenitors
scv.tl.score_genes_cell_cycle(adata)
scv.pl.scatter(adata, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95])

s_genes, g2m_genes = scv.utils.get_phase_marker_genes(adata)
s_genes = scv.get_df(adata[:, s_genes], 'spearmans_score', sort_values=True).index
g2m_genes = scv.get_df(adata[:, g2m_genes], 'spearmans_score', sort_values=True).index
kwargs = dict(frameon=False, ylabel='cell cycle genes')
scv.pl.scatter(adata, list(s_genes[:2]) + list(g2m_genes[:3]), **kwargs)
scv.pl.velocity(adata, ['HELLS', 'TOP2A'], ncols=2, add_outline=True)
```

经验教训：生信分析所用的很多软件由于维护者不会及时更新，对包与包之间的依赖及版本要求很高，一些新的python及R的版本不一定能满足大部分包，尤其是一些上古神包的运行要求，在安装加载各个过程中都可能出错。有大神推荐在python3.8版本下安装单细胞分析软件会比较稳定，总结就是装旧不装新（软件包），不随便更新（python或R语言及包）。

最后尤其感谢**Jimmy**老师等开源代码及生信教程维护者，让我从一个生信小白到能进行一些简单的个性化分析，也希望通过这样的总结帮助更多人，上述代码也会放在Github上：https://github.com/PrinceWang2018/scvelo_py3.10