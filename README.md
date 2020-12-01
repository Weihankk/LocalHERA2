## 注意
**LocalHERA2将BWA替换成了minimap2， 速度有了质的提升，但与BWA比对结果的差异还未知，目前BWA还在运行中，根据HERA论文附件说明是完全可以用minimap的。**

**更新日期：2020.12.1**

# LocalHERA2
This is a new respository for running HERA in local environment

![git](https://img.shields.io/badge/HERA-Local-brightgreen) 

Local Highly Efficient Repeat Assembly 2 

## 项目介绍
将HERA流程拆分，并使用R重写了本地运行流程。

HERA项目源地址：https://github.com/liangclab/HERA

LocalHERA旧版地址（**已弃用**）：https://github.com/Weihankk/LocalHERA

如果使用本流程请引用HERA原文：

Du, H., & Liang, C. (2019). Assembly of chromosome-scale contigs by efficiently resolving repetitive sequences with long reads. Nature communications, 10(1), 1-10.

本人邮箱：whzhang@webmail.hzau.edu.cn

本人QQ：97578011


## 使用方法
### 所需数据
- CANU/FALCON等软件组装的contigs.fasta文件
- CANU/FALCON等软件使用三代自校正的corrected.fasta文件

### 准备
```
mkdir RunHERA2
cd RunHERA2
git clone https://github.com/Weihankk/LocalHERA2.git
cd LocalHERA2
vi LocalHERA2_minimap.R  # 修改HERA运行脚本中的软件路径，先将就着用，以后会整合到参数文件中 ...
# 修改第94行、125行的minimap路径，以及251-261行的dalign路径为自己服务器上的路径即可
cd ..
```
Now you are in `/RunHERA2` directory

### 配置参数文件
- 首先创建一个HERA工作目录
```
mkdir HERA2_TEST
cp ../LocalHERA2/LocalHERA2_Parameters.R ./
```
- 修改参数文件`LocalHERA2_Parameters.R`, 主要是工作目录，contig fasta， correct fasta，其他参数没搞懂我也就没变动
```
vi LocalHERA2_Parameters.R
```

### Run
```
Rscript LocalHERA2_minimap.R LocalHERA2_Parameters.R
```

### 加速
- `LocalHERA2_Parameters.R` 中`cl <- makeCluster(5)` 表示同时运行5个dalign任务，可以根据自己服务器情况增加，但不得超过单节点的总线程数。
-  94行和125行minimap2命令中的线程数也可以自行改变


