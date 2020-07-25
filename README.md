# LocalHERA2
This is a new respository for running HERA in local environment

![git](https://img.shields.io/badge/HERA-Local-brightgreen) 

Local Highly Efficient Repeat Assembly 2 

## 项目介绍
将HERA流程拆分，并使用R重写了本地运行流程脚本，**修改了之前的一些小问题**。

HERA项目源地址：https://github.com/liangclab/HERA

LocalHERA源地址（**已弃用**）：https://github.com/Weihankk/LocalHERA

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
mkdir RunHERA
cd RunHERA
git clone https://github.com/liangclab/HERA.git
git clone https://github.com/Weihankk/LocalHERA2.git
cd LocalHERA2
cp ../HERA/HERA.zip ./
unzip HERA.zip
cd HERA
chmod 755 *
cp ../../HERA/*fasta ./
cd ../
```
Now you are in `/RunHERA` directory

**HERA部署成功**

### 配置参数文件
- 首先创建一个HERA工作目录
```
mkdir HERA_TEST
```
- 修改参数文件`LocalHERA_Parameters.R`
```
vi LocalHERA_Parameters.R
```
- 按照注释内容分别将各种数据和工具指定好，此处相比HERA原版只增加了前4行：`set working directory`、`set bwa-0.7.10`

### Step 1




