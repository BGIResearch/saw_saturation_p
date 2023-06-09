# saturation

Calculate sequencing saturation in SAW

## 程序版本

本文档对应程序版本为 **v2.3.5**

## 运行环境要求

linux系统环境

python3要求安装,且在环境变量PATH中,需要安装python库包括:

* numpy
* pandas
* matplotlib
* scipy

## 输入文件要求

* 原始基因表达文件. 名称为: raw_barcode_gene_exp.txt **兼容格式有表头/无表头,表头需按特定格式,五列数据以空格分隔**
* 组织区域坐标文件. 名称为: SN.gem.gz或SN.gef **前者要求格式有注释行,有表头,四列数据以tab分隔**
* 比对统计文件. 名称为: SN.barcodeMap.stat **列以tab分隔,包含以 total_reads: 开头的行**
* 注释统计文件. 名称为: SN.summary.stat **要求有表头,列以tab分隔**

## 输出文件说明

* 测序饱和度统计文件.
  文件共九列,以tab分隔;第一行为表头,示例如下:
  | 抽样比例 | bin1 total reads | bin1 饱和度值 | bin1 基因数中值 | bin1 uniq reads | bin200 total reads | bin200 饱和度值 | bin200 基因数中值 | bin200 uniq reads |
  | -------- | ---------------- | ------------- | --------------- | --------------- | ------------------ | --------------- | ----------------- | -------------- |
  | #sample  | bar_x            | bar_y1        | bar_y2          | bar_umi         | bin_x              | bin_y1          | bin_y2            | bin_umi        |
  | 0.05     | 107416           | 0.5978346     | 1               | 43199           | 107416             | 0.5981511       | 13                | 51             |
  | 0.1      | 214832           | 0.718608      | 1               | 60452           | 214832             | 0.7188687       | 19                | 69             |
  | 0.2      | 429665           | 0.8042382     | 1               | 84112           | 429665             | 0.8044639       | 29                | 91             |
  | 0.3      | 644498           | 0.8399685     | 1               | 103140          | 644498             | 0.8402043       | 37                | 109            |
  | 0.4      | 859330           | 0.8606996     | 1               | 119705          | 859330             | 0.8609568       | 43                | 124            |
  | 0.5      | 1074163          | 0.8747304     | 1               | 134560          | 1074163            | 0.8749957       | 49                | 137            |
  | 0.6      | 1288996          | 0.8849469     | 1               | 148303          | 1288996            | 0.8852138       | 55                | 148            |
  | 0.7      | 1503828          | 0.8928634     | 1               | 161115          | 1503828            | 0.8931374       | 58                | 160            |
  | 0.8      | 1718661          | 0.8992623     | 1               | 173134          | 1718661            | 0.8995398       | 62                | 170            |
  | 0.9      | 1933494          | 0.9044445     | 1               | 184756          | 1933494            | 0.9047227       | 64                | 180            |
  | 1        | 2148327          | 0.9088984     | 1               | 195716          | 2148327            | 0.9091828       | 66                | 189            |
* 饱和度图片. 所有图片的x轴都是相同的,读取统计文件第二列(bin1 total reads),并按比例还原回total reads,即最后一个采样点的x轴坐标对应total reads(即可视化页面的total reads)
  * bin1图片. 共两幅图,分别以第三列(bin1 饱和度值),第四列(bin1 基因数中值)作为y轴
  * bin200图. 共三幅图,前两幅分别以第七列(bin200 饱和度值),第八列(bin200 基因数中值)作为y轴; 第三幅图以最后一列(bin200 uniq reads)作为y轴,并对该统计结果做对数拟合,将拟合公式输出到统计图上,如果最大取样量下的值小于设定的阈值5000,并且拟合曲线的R2 大于等于0.9,就画出拟合曲线的延长线,一直延长到阈值处,并在图上画出阈值处的坐标值

## 程序处理逻辑

饱和度值计算公式为 *1-(uniq reads/total reads)*

* 准备. 
  * 读取比对统计文件,获取total reads
  * 读取注释统计文件,获取anno reads
  * 计算比例anno reads/total reads,用于画图时还原x轴最大坐标为total reads
* 读tissue数据文件. 
  * 将所有unique 坐标提取出来,注意这里的坐标需要加上offset
  * 将坐标取bin,得到bin200的unique 坐标
  * 随机采样5%的bin200 unique坐标,并还原回bin1的坐标,用于过滤表达矩阵中的数据
* 读表达矩阵文件.
  * 使用tissue下的采样后的bin1坐标过滤数据,构建(x,y,gene,mid)的列表
  * 同时累加所有count值,得到anno reads
* 饱和度计算.
  * 将前一列表顺序打乱
  * 按照{0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}的采样间隔,顺序处理数据
  * 对每一个采样点,分别计算bin1/bin200下的total reads/饱和度值/基因数中值等,输出统计结果文件
    * **饱和度值**. 计算所有bin下的uniq reads(即gene,mid都为uniq)和total reads,使用公式1-(uniq reads/total reads)计算得出; bin1 bin200算法相同
    * **基因数中值**. 计算每个bin下的uniq gene个数,取中值; bin1 bin200算法相同
    * **uniq reads**. bin1的uniq reads为所有bin1下的uniq reads(即gene,mid都为uniq); bin200的uniq reads为所有bin200下的uniq reads(bin1的x y坐标,gene和mid都为uniq)
* 绘图. 分别绘制bin1/bin200的统计图
