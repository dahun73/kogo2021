
# ***ATTENTION***
아래 코드는 <br>
2021년 07월 07-08 일 통계유전체학회에서 실시하는<br>
인체마이크로바이옴분석 실습 수업을 위해 준비한 것입니다. <br>
본 메뉴얼의 저작권은 [QIIME2](https://docs.qiime2.org/2021.4/citation/) 에 있습니다. <br>
강의자료 제작: <span style="color:blue">강북삼성병원 **김한나**</span> (문의: 147942@hanmail.net) 
----


파란 글씨는 클릭이 가능하며, 오른쪽 클릭으로 새창으로 열기로 여시면, 해당 웹페이지로 연결됩니다.



# "Atacama soil microbiome" tutorial
- [튜토리얼 링크](https://docs.qiime2.org/2021.4/tutorials/atacama-soils/)를 참고하시되, input 파일 포멧이 튜토리얼과 차이가 있습니다.
- 튜토리얼에는 EMP 데이터셋(multiplexed paired-end)으로 Importing 후 Demultiplexing => Denoising 이 진행되지만,
- 본 실습에서는 일반적으로 연구자가 처음 받게되는 demultiplexed paired-end 데이터셋(Casava 1.8 demultiplexed paired-end sequences)으로 Importing 후 바로 Denoising 단계로 진행합니다.

- 튜토리얼은 두가지 목적으로 만들어졌습니다.

---
1. Paired-end read analysis
    - Importing
    - Demultiplexing (생략. 이미 demultiplexed 된 데이터 이용)
    - Denoising
    - 결과 : Feature table & Feature sequences

2. Self-guided exercise
    - QIIME2 에 좀 더 익숙해지기 위하여, single-end read analysis 인 ["Moving pictures" tutorial](https://docs.qiime2.org/2021.4/tutorials/moving-pictures/) 복습 해보기
    - [튜토리얼](https://docs.qiime2.org/2019.4/tutorials/atacama-soils/) 마지막 부분의 questions 에 대하여 답해보기
---

> #### 실습 파일 (이메일 링크 다운로드)
- <span style="color:red"> 바탕화면의 KOGO_2021/Day2/atacama_tutorial </span> 폴더내
    - Input files
    	- **casava-18-paired-end-demultiplexed 폴더** : 61 samples 의 paired-end fastq 파일 (V4 region)
        - **sample-metadata_ata.tsv** : phenotype 파일
    - reference DB (Naive Bayes classifiers로 trained 완료)
        - Silva version 138, 99% OTUs, V4 region : **silva-138-99-515-806-nb-classifier.qza**
            
- 아래 실습 단계에서 생성되는 **1) ~ 21)** 의 모든 결과는 **atacame_answer 폴더** 에 있으며, 각자 실행하여 얻은 결과와 비교 가능합니다.
  	  	
   
---


## Contents ##


* [QIIME 2 시작하기](#qiime-2-시작하기)
    * [1. Importing](#1-importing)
    * [2. Denoising using DADA2 and Feature Table 만들기](#2-denoising-using-dada2-and-feature-table-만들기)
        * [DADA2 실행](#dada2-실행)
        * [QC Summary 파일 Visualization](#qc-summary-파일-visualization)
    * [3. Feature Table and Feature Data Summaries](#3-feature-table-and-feature-data-summaries)
        * [Feature Table Summary 파일 만들기](#feature-table-summary-파일-만들기)
        * [Representative Sequences 확인](#representative-sequences-확인)
     * [4. Filtering Data](#4-filtering-data)
     * [5. Phylogenetic Diversity 분석을 위한 Phylogenetic Tree 만들기](#5-phylogenetic-diversity-분석을-위한-phylogenetic-tree-만들기)
     * [6. Alpha and Beta Diversity Analysis](#6-alpha-and-beta-diversity-analysis)
         * [Core Analysis](#core-analysis)
         * [Alpha Diversity 그룹간 비교](#alpha-diversity-그룹간-비교)
         * [Alpha Diversity 연속변수 분석](#alpha-diversity-연속변수-분석)
         * [Beta Diversity 그룹간 비교](#beta-diversity-그룹간-비교)
     * [7. Taxonomic Analysis](#7-taxonomic-analysis)
         * [Classification](#classification)
         * [Taxonomic Composition 에 대해 Bar Plot 그리기](#taxonomic-composition-에-대해-bar-plot-그리기)
     * [8. Group 간 Differential Abundance Test](#8-group-간-differential-abundance-test)
         * [Taxa Collapsing (Taxa level 별로 분석하기)](#taxa-collapsing-taxa-level-별로-분석하기)
          * [Phylum 수준에서 분석 (L2)](#phylum-수준에서-분석-l2)
         * [Genus 수준에서 분석 (L6)](#genus-수준에서-분석-l6)
     * [9. Exporting](#9-exporting)
    
        
   
# QIIME 2 시작하기

- Mac 또는 Window 에서 qiime2 를 활성화 시킵니다.

- 실습폴더(바탕화면)로 이동합니다. (사용자마다 실습폴더 위치는 다를 수 있습니다)
    - 예) hanna 컴퓨터의 바탕화면에 있는 실습폴더 

```sh
# 예) Mac, hanna 컴퓨터의 바탕화면에 있는 실습폴더로 이동
cd /Users/hanna/Desktop/KOGO_2021/Day2/atacama_tutorial

# 예) Windows, Ubuntu 에서 바탕화면에 있는 실습폴더로 이동
cd /mnt/c/Users/hanna/Desktop/KOGO_2021/Day2/atacama_tutorial

# 제대로 들어왔는지 현재폴더위치를 확인합니다.
pwd

# QIIME2 활성화 시키기
conda activate qiime2-2021.4
```
</br></br>
# 1. Importing
### Casava 1.8 paired-end demultiplexed fastq 포멧

- [Casava 1.8 demultiplexed](http://illumina.bioinfo.ucr.edu/ht/documentation/data-analysis-docs/CASAVA-FASTQ.pdf/view) 포멧은 아래 5가지 이름이 underscore로 연결된 파일명 형태이며, Illumina Miseq 으로 paired-end 시퀀싱 시 일반적으로 연구자가 받게되는 파일명 형태입니다.
    - 예) L2S357_15_L001_R1_001.fastq.gz, L2S357_15_L001_R2_001.fastq.gz
    1. sample identifier
    2. barcode sequence or a barcode identifier
    3. lane number
    4. direction of the read (i.e. R1 or R2)
    5. set number
    
- 실험 의뢰기관으로부터 받은 fastq 파일 이름이, 위와 같은 5개 구분형태가 아니라면, 되도록 위와같은 형태로 요구하여 받으시는게,
튜토리얼을 따라가며 작업하시는데 도움이 됩니다.


---
> #### 참고

- 만약, Casava 1.8 demultiplexed 형태로 받을 수 없다면, 파일명을 Casava 포멧으로 전체 변경하거나, "Fastq manifest" 포멧을 이용하여 import 해야합니다 ([관련설명: Fastq manifet format](https://docs.qiime2.org/2021.4/tutorials/importing/))


---



### **Importing** 하기
pwd 했을때, 현재 실습폴더 (Desktop/KOGO_2021/Day2/atacama_tutorial)에 위치해 있다는 가정하에 아래 명령어를 그대로 copy&paste 하시기 바랍니다.
- Casava 1.8 포멧으로 네이밍된 fastq 파일들 모두(모든 샘플) 한 폴더에 넣어두고, 아래 input-path 에 해당 폴더명을 입력합니다.
```sh
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path casava-18-paired-end-demultiplexed \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path casava_pe_demux.qza
```
- output 파일: **01) casava_pe_demux.qza** 생성확인
- casava_pe_demux.qza 파일을 Visualtizations
```sh
qiime demux summarize \
--i-data casava_pe_demux.qza \
--o-visualization casava_pe_demux.qzv
```
- output 파일: **2) casava_pe_demux.qzv** 생성확인
- 위 .qzv 파일을 [QIIME2view](https://view.qiime2.org) 에서 열기. 
   - Overview : 샘플 데이터의 총 시퀀스(리드) 통계 및 샘플 별 시퀀스 수
   - Interactive Quality Plot : Forward 와 Reverse Reads 의 Quality scores 를 Plot 에서 확인하고, QC 단계의 조건들을 고민함
- fastq 파일 이외에 포멧을 QIIME2 에서 이용하고 싶다면, 

</br></br>


# 2. Denoising using DADA2 and Feature Table 만들기

DADA2 에서는, 시퀀스의 QC 뿐만 아니라, paire-end 시퀀싱 Merging 및 Chimera 제거도 동시에 실행됩니다.
- casava_pe_demux.qzv 파일을 통해 직접 눈으로 보고 확인한 read quality 를 기반으로 trimming과 truncation 결정
    - --p-trim-left-f 과 --p-trim-left-r
    - --p-trunc-len-f 과 --p-trunc-len-r
- truncation 시, paired-end 의 overlap 길이를 반드시 고려하여, 충분한 overlap이 유지되도록 주의요함
- 실습파일: 2x150bp paired end 리드로서, V4 영역(약 250bp) 디자인. read triming 전혀 하지 않아도 50bp 오버랩. 최소 30bp 오버랩구간을 남기려면, 현재 20bp 트리밍이 최대 가능.

> #### DADA2 실행
```sh
qiime dada2 denoise-paired \
--i-demultiplexed-seqs casava_pe_demux.qza \
--p-trim-left-f 10 \
--p-trim-left-r 10 \
--p-trunc-len-f 150 \
--p-trunc-len-r 150 \
--o-table table-dada2.qza \
--o-representative-sequences rep-seqs.qza \
--o-denoising-stats stats-dada2.qza
```

- **03) table-dada2.qza** 생성확인
- **04) rep-seqs.qza** 생성확인
- **05) stats-dada2.qza** 생성확인

**- DADA2로 생성된 파일 중, table-dada2.qza 파일은, Feature Table 로서, Downstream analysis 에 계속 사용될 core input file 이니 기억하세요!!**

> #### QC Summary 파일 Visualization 
```sh
qiime metadata tabulate \
--m-input-file stats-dada2.qza \
--o-visualization stats-dada2.qzv
```
- **06) stats-dada2.qzv** 생성확인
- 위 .qzv 파일을 [QIIME2view](https://view.qiime2.org) 에서 열기. 
- 최초 raw data 가 DADA2 작업 후, 최종 어느정도 filtering되었는지 확인가능
    - 샘플별 최초 시퀀스 리드수(input), QC로 filter된 후(filtered, denoised), paired-end 시퀀스 merging후(merged), chimera 제거후(non-chimeric)의 리드수 확인가능
    - filtered reads 수가 너무 적게 남았거나하면, trimming 이나 truncation parameter 수를 조정하여 DADA2 를 재실행함
</br></br>
# 3. Feature Table and Feature Data Summaries

> #### Feature Table Summary 파일 만들기

```sh
qiime feature-table summarize \
--i-table table-dada2.qza \
--o-visualization table-dada2.qzv \
--m-sample-metadata-file sample-metadata_ata.tsv
```
- **07) table-dada2.qzv** 생성확인
- 위 .qzv 파일을 [QIIME2view](https://view.qiime2.org) 에서 열기. 
   
   - Overview
        - Sample수, features 수, 총 리드수(frequency)
        - 샘플별 frequency(리드수)
        - Feture별 frequency
    - Interactive Sample Detail
        - phenotype 데이터 정보에 기반하여, 샘플별 리드수, 변수(그룹)별 리드수 확인가능
        - Diversity 분석시 rarefaction 위한 sampling depth 기준 설정을, 여기서 확인함
    - Feature Detail
        - Amplicon Sequence Variants (ASV)당 frequency

> #### Representative Sequences 확인
```sh
qiime feature-table tabulate-seqs \
--i-data rep-seqs.qza \
--o-visualization rep-seqs.qzv
```
- **08) rep-seqs.qzv** 생성확인
- 위 .qzv 파일을 [QIIME2view](https://view.qiime2.org) 에서 열기. 
   
   - Sequence Length 통계 : ASVs 마다의 평균 시퀀스 length
    - Sequence Lenghs의 평균 분포
    - Sequence Table : Features의 실제 시퀀스

</br></br>
# 4. Filtering Data

위 생성된 03) tabble-dada2.qzv 파일을 viewer 에서 확인 후, 리드 수가 현저히 낮은 샘플을 제거하기로 함

```sh
qiime feature-table filter-samples \
--i-table table-dada2.qza \
--p-min-frequency 100 \
--o-filtered-table filtered_100_table.qza
```
- **09) filtered_100_demux.qza** 파일 생성 확인
- 이 filtered_100_demux.qza 가 무엇이 바뀌었는지 **10) filtered_100_demux.qzv** 로 만들어서 뷰어로 확인하기 (샘플 수, 리드 수)

```sh
qiime feature-table summarize \
--i-table filtered_100_table.qza \
--o-visualization filtered_100_table.qzv
--m-sample-metadata-file sample-metadata_ata.tsv
```
- sample 또는 feature 필터 기능은, 실제 많이 사용되는 기능이니 필요 시, [QIIME2홈페이지-Filtering Data](https://docs.qiime2.org/2021.4/tutorials/filtering/) 에서 필요한 부분을 응용하시기 바랍니다.

</br></br>

# 5. Phylogenetic Diversity 분석을 위한 Phylogenetic Tree 만들기

```sh
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs.qza \
--o-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza
```
- **11) aligned-rep-seqs.qza** 생성확인
- **12) masked-aligned-rep-seqs.qza** 생성확인
- **13) unrooted-tree.qza** 생성확인
- **14) rooted-tree.qza** 생성확인 => 이 rooted-tree.qza 파일이 phylogenetic diversity 분석에 사용됩니다.

</br></br>

# 6. Alpha and Beta Diversity Analysis

Diversity 분석은, QIIME2의 "diversity" 라는 plugin 을 사용하며, core-metrics-phylogenetic 방법으로, alpha 와 beta diversity metrics를 한번에 생성가능합니다.

- Alpha diversity
    - Shannon’s diversity index (a quantitative measure of community richness)
    - Observed Features (a qualitative measure of community richness)
    - Faith’s Phylogenetic Diversity (a qualitiative measure of community richness that incorporates phylogenetic relationships between the features)
    - Evenness (or Pielou’s Evenness; a measure of community evenness)
- Beta diversity
    - Jaccard distance (a qualitative measure of community dissimilarity)
    - Bray-Curtis distance (a quantitative measure of community dissimilarity)
    - unweighted UniFrac distance (a qualitative measure of community dissimilarity that incorporates phylogenetic relationships between the features)
    - weighted UniFrac distance (a quantitative measure of community dissimilarity that incorporates phylogenetic relationships between the features)

- core-metrics-phylogenetic 방법에서는, FeatureTable[Frequency] 즉, 위에서 생성한 **09) filtered_100_table.qza** 파일을 user-specified depth 로 rarefying 합니다.
- 이를 위하여, **10) filtered_100_demux.qzv** 을 [QIIME2view](https://view.qiime2.org) 의 "Interactive Sample Detail"에서 확인한대로, 적절한 depth를 정합니다. 정한 depth 로 replacement 없이 subsampling을 하며, 기준 depth 이하의 샘플은 분석에서 제외됩니다.


> #### Core Analysis
- --p-sampling-depth 명령어와 함께 해당 depth 숫자를 반드시 적습니다.
```sh
qiime diversity core-metrics-phylogenetic \
--i-phylogeny rooted-tree.qza \
--i-table filtered_100_table.qza \
--p-sampling-depth 142 \
--m-metadata-file sample-metadata_ata.tsv \
--output-dir core-metrics-results
```
- 현재 실습폴더에 **15) core-metrics-results** 폴더가 새로 생성된 것을 확인하세요.
- 위 명령어(--output-dir) 실행 전에, 동일이름의 폴더가 존재하면 명령이 실행되지않고 에러가 남. 반드시 작업폴더에 존재하지 않는 New Name 을 --output-dir 뒤에 적습니다.
- core-metrics-results 폴더로 들어가보세요.
```sh
cd core-metrics-results
```
- 폴더에 어떤 파일들이 생성되었는지 확인해보세요.
```sh
ls
```
- 아래와 같이 확인되며, alpha diversity indices 와 beta diversityi indices 가 섞여서 나열되어있으니, 위 섹션에서 alpha 와 beta 에 각각 해당되는 것을 확인하고 다음 분석에 이용하세요.
    - bray_curtis_distance_matrix.qza
    - bray_curtis_emperor.qzv
    - bray_curtis_pcoa_results.qza
    - evenness_vector.qza			
    - faith_pd_vector.qza			
    - jaccard_distance_matrix.qza		
    - jaccard_emperor.qzv
    - jaccard_pcoa_results.qza
    - observed_features_vector.qza
    - rarefied_table.qza
    - shannon_vector.qza
    - unweighted_unifrac_distance_matrix.qza
    - unweighted_unifrac_emperor.qzv
    - unweighted_unifrac_pcoa_results.qza
    - weighted_unifrac_distance_matrix.qza
    - weighted_unifrac_emperor.qzv
    - weighted_unifrac_pcoa_results.qza
- core 분석으로 생성된 파일들은 연구목적의 그룹간 diversity 비교분석(통계)을 위한 input files 이라 생각하시면 되며, 위의 .qza 파일들 자체만으로는 아무것도 알 수 없습니다.
- 단, core 분석의 beta-diversity 관련 결과 중, *_emperor.qzv 파일은, [QIIME2view](https://view.qiime2.org) 에서 바로 확인가능합니다.

    - bray_curtis_emperor.qzv를 [QIIME2view](https://view.qiime2.org) 에서 확인해보세요. 

- core-metrics-results 폴더에서 다시 실습폴더(한단계 상위폴더)로 돌아옵니다.
```sh
# 한단계 상위폴더로 올라가기
cd ..

# 현재 위치가 바탕화면의 QIIME2실습폴더 아래 atacama_tutorial 폴더 위치가 맞는지 확인하기
pwd

```
 
> #### Alpha Diversity 그룹간 비교

- input file 의 경로와 output file 의 경로에 주의하여 아래 명령어를 실행합니다.
- Faith's PD (phylogenetic diversity)
```sh
qiime diversity alpha-group-significance \
--i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
--m-metadata-file sample-metadata_ata.tsv \
--o-visualization core-metrics-results/faith-pd-group-significance.qzv
```
- **15) core-metrics-results** 에 faith-pd-group-significance.qzv 신규 생성 확인
- 해당 .qzv 파일을 [QIIME2view](https://view.qiime2.org) 에서 열어보세요. 
     
- Alpha Diversity Boxplots 확인
    - Column 탭을 조절하여 비교하고싶은 변수(그룹)별 비교가능
    - Kruska-Wallis 통계량 확인
    - 모든 plots 은 SVG 형식으로 다운 가능하며, TSV나 CSV 다운로드하여 직접 plot을 새로 그릴 수도 있음

- Pielou's Evenness
```sh
qiime diversity alpha-group-significance \
--i-alpha-diversity core-metrics-results/evenness_vector.qza \
--m-metadata-file sample-metadata_ata.tsv \
--o-visualization core-metrics-results/evenness-group-significance.qzv
```
- **15) core-metrics-results 폴더** 에 evenness-group-significance.qzv 신규 생성 확인
- 해당 .qzv 파일을 [QIIME2view](https://view.qiime2.org) 에서 열어보세요. 
- Alpha diversity 의 다른 index 도 돌려보세요
    - observed_otus_vector.qza
    - shannon_vector.qza


> #### Alpha Diversity 연속변수 분석
- 그룹(categorical)간 비교 뿐만 아니라, 나이 등의 연속변수(numeric)도 분석가능합니다.
- sample metadata (phenotype)의 2행에, 변수별 성격 (categorical/numeric)가 미리 입력되어 있어야하며, 연속변수 분석은 특정 하나의 변수 지정이 아닌, metadata에 포함된 모든 연속변수에 대하여 진행됩니다.
- --p-method 는 spearman 방법과 pearson 둘 중 하나를 사용할 수 있습니다.

```sh
qiime diversity alpha-correlation \
--i-alpha-diversity core-metrics-results/shannon_vector.qza \
--p-method spearman \
--m-metadata-file sample-metadata_ata.tsv \
--o-visualization core-metrics-results/qt_shannon.qzv
```
- **15) core-metrics-results 폴더** 에 qt_shannon.qzv 신규 생성 확인
- 해당 .qzv 파일을 [QIIME2view](https://view.qiime2.org)에서 열어보세요. 
- 결과보기
    - Alpha Correlation : Column 에서 해당 연속변수를 선택하면 각각의 연속변수에 대한 correlation plot 확인가능
    - spearman (또는 pearson)통계 결과 확인가능
    - plot은 SVG 파일로 다운가능하며, 통계분석의 input 파일 역시 TSV 파일로 다운로드 가능함


> #### Beta Diversity 그룹간 비교 
- Unweighted UniFrac distance (phylogenetic binary 분석)
```sh
qiime diversity beta-group-significance \
--i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
--m-metadata-file sample-metadata_ata.tsv \
--m-metadata-column Vegetation \
--o-visualization core-metrics-results/unweighted-unifrac-vegetation-significance.qzv \
--p-pairwise
```
- **15) core-metrics-results** 에 unweighted-unifrac-vegetation-significance.qzv 신규 생성 확인
- 해당 .qzv 파일을 [QIIME2view](https://view.qiime2.org) 에서 열어보세요. 
- 결과 보기
    - Overview : PERMAVONA 통계량 확인
    - Group significance plots : Distance to 기준그룹
    - Pairwise permanova results : 세그룹 이상 시, pairwise 통계분석결과
    - 모든 plots 은 PDF 다운 가능하며, TSV나 CSV 다운로드하여 직접 plot을 새로 그릴 수도 있음

- Weighted UniFrac distance (phylogenetic abundance 분석) 
```sh
qiime diversity beta-group-significance \
--i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
--m-metadata-file sample-metadata_ata.tsv \
--m-metadata-column Vegetation \
--o-visualization core-metrics-results/weighted-unifrac-vegetation-significance.qzv \
--p-pairwise
```
- **15) core-metrics-results** 에 weighted-unifrac-vegetation-significance.qzv 신규 생성 확인
- 해당 .qzv 파일을 [QIIME2view](https://view.qiime2.org) 에서 열어보세요. 
- Beta diversity 의 다른 index 도 돌려보세요
    - bray_curtis_distance_matrix.qza (non-phylogenetic/abundance)
    - jaccard_distance_matrix.qza (non-phylogenetic/binary) 
 
 
- PCoA plot 확인
    - unweighted_unifrac_emperor.qzv 을 [QIIME2view](https://view.qiime2.org) 에서 열어보세요.*_emperor.qzv 파일들은 core_analysis 분석 시 자동으로 qzv 파일이 생성되어 있습니다.
 
</br></br>
# 7. Taxonomic Analysis

> #### Classification  
- Taxonomic 분석을 위해서는, 기존 만들어놓은 feature table ( **09) filtered_100_table.qza** )과 representative sequence ( **08) rep-seqs.qza**) 파일 외에, taxonomy (features 의 이름) 정보가 필요합니다.
- 그러므로 taxonomic 분석 전에, **08) rep-seqs.qza** 파일에 Machine Learning 방법으로 train된 reference  DB (QIIME2 홈페이지에서 [다운로드](https://docs.qiime2.org/2021.4/data-resources/) 가능)의 taxonomy 정보를 붙이는 작업을 아래와 같이 진행합니다.
    - Silva DB 를 이용한 trained data : silva-138-99-515-806-nb-classifier.qza
    - Classification 작업은 시간이 조금 오래 걸립니다.
```sh
qiime feature-classifier classify-sklearn \
--i-classifier silva-138-99-515-806-nb-classifier.qza \
--i-reads rep-seqs.qza \
--o-classification silva138_99_taxonomy.qza
```
- 현재 실습폴더에 **16) silva138_99_taxonomy.qza** 폴더가 새로 생성된 것을 확인하세요.

- representative sequences 데이터에 taxonomy 가 붙은 것을 확인하기 위해, 위 생성된 .qza 파일을 Visualization
```sh
qiime metadata tabulate \
--m-input-file silva138_99_taxonomy.qza \
--o-visualization silva138_99_taxonomy.qzv
```
- 현재 실습폴더에 **17) silva138_99_taxonomy.qzv** 파일이 새로 생성된 것을 확인하세요.
- 해당 .qzv 파일을 [QIIME2view](https://view.qiime2.org) 에서 열어보세요.


> #### Taxonomic Composition 에 대해 Bar Plot 그리기 
```sh
qiime taxa barplot \
--i-table filtered_100_table.qza \
--i-taxonomy silva138_99_taxonomy.qza \
--m-metadata-file sample-metadata_ata.tsv \
--o-visualization taxa-bar-plots_silva.qzv
```
- 현재 실습폴더에 **18) taxa-bar-plots_silva.qzv** 파일이 새로 생성된 것을 확인하세요.
- 해당 .qzv 파일을 [QIIME2view](https://view.qiime2.org) 에서 열어보세요. 

- Taxonomic Level 탭에서 변경가능 : L1 (Kindom), L2 (Phylum), L3 (Class), L4 (Order), L5 (Family), L6 (Genus), L7 (Species)
- 변수(그룹)에 따라 표시 가능
- Plot 은 SVG 파일로 다운가능하며, CSV 파일로도 다운로드 가능하므로, R 등의 다른 프로그램으로 ploting 가능합니다.

</br></br>
# 8. Group 간 Differential Abundance Test 
- Differential abundance testing 을 위한 QIIME2 plugins 은 아래와 같이 3종류가 있으며, 튜토리얼에서는 composition plugin 인 **ANCOM**을 이용함
    - ANCOM ([논문](https://pubmed.ncbi.nlm.nih.gov/26028277/))
    - gneiss ([설명](https://docs.qiime2.org/2021.4/tutorials/gneiss/?highlight=gneiss))
    - ALDEx2 ([plugin 추가설치 필요](https://library.qiime2.org/plugins/q2-aldex2/24/))
- ANCOM은 features 중에 약 25% 미만에서 차이난다는 가정하에 분석하는 것이므로, 만약 내 데이터의 그룹간에 더 많은 features 가 차이날꺼라 생각한다면 ANCOM을 사용해서는 안됨 (Type I/ Type II 에러 모두 증가)
    - 예) Gut vs. Skin, Gut vs. Oral
- 일반적인 분석에서 gut 에서만 또는 skin 에서만 연구를 진행하는 경우가 많고, 같은 부위에서는 diversity 분석에서도 극명한 차이를 나타내지 않으므로, ANCOM 을 tutorial 에서 진행함

- Taxonomic Analysis 분석을 위해서는, feature table ( **09) filtered_100_table.qza** )와 metadata (**sample-metadata_ata.tsv**) 가 필요합니다.

    
> #### Taxa Collapsing (Taxa level 별로 분석하기)
- L2-L7 까지 collapsing 가능
    - 예) L2 Phylume level 로 collapsing 하기
- 다른 taxa level을 만들때도 아래와 같이 3단계(collapse, add-pseudocount, ancom)를 동일하게 진행해야합니다.
- collapse 단계 없이, filtered_100_table.qza 로 바로 add-pseudocount 부터 들어가면 ASV level 로 분석하게 되며, 분석 시간이 매우 오래 걸릴 수 있습니다.

#### Phylum 수준에서 분석 (L2)
```sh
qiime taxa collapse \
--i-table filtered_100_table.qza \
--i-taxonomy silva138_99_taxonomy.qza \
--p-level 2 \
--o-collapsed-table table_L2.qza
```
- 현재 실습폴더에 **19) table_L2.qza** 파일이 새로 생성된 것을 확인하세요.

```sh
qiime composition add-pseudocount \
--i-table table_L2.qza \
--o-composition-table comp_table_L2.qza
```
- 현재 실습폴더에 **20) comp_table_L2.qza** 폴더가 새로 생성된 것을 확인하세요.

```sh
qiime composition ancom \
--i-table comp_table_L2.qza \
--m-metadata-file sample-metadata_ata.tsv \
--m-metadata-column Vegetation \
--o-visualization ancom-vegetation_L2.qzv
```
- 현재 실습폴더에 **21) ancom-vegetation_L2.qzv** 파일이 새로 생성된 것을 확인하세요.
- 해당 .qzv 파일을 [QIIME2view](https://view.qiime2.org) 에서 열어보세요. 

  - 현재 분석에서는 total 16 개의 phyla 중, 가티 8개의 phyla와 다르게 두그룹사이에서 유의미하게 차이난 2개의 phyla 가 W값 8을 가지며 TRUE 결과를 제시하고 있습니다.
  - Volcano plot 과 Percentile abundances 결과로부터 "yes" 그룹이 "no" 그룹보다 많이 가지는 균임을 확인할 수 있습니다.
    - ANCOM Volcano Plot : X 축의 0을 기준으로 우측은 비교군에서 증가, 좌측은 감소. Y 축은 W값을 의미하며, taxa level 에 따라 상대적인 값임
    - ANCOM statistical results : ANCOM 분석 결과 유의한 taxa 만 결과로 나타남. 유의미한게 없을 시, "No significant features found"
  - (유의미한 결과 포함한) 전체 결과는 "Download table as TSV" 클릭하여 다운로드 가능


#### Genus 수준에서 분석 (L6)
```sh
qiime taxa collapse \
--i-table filtered_100_table.qza \
--i-taxonomy silva138_99_taxonomy.qza \
--p-level 6 \
--o-collapsed-table table_L6.qza
```
- 현재 실습폴더에 **22) table_L6.qza** 폴더가 새로 생성된 것을 확인하세요.

```sh
qiime composition add-pseudocount \
--i-table table_L6.qza \
--o-composition-table comp_table_L6.qza
```
- 현재 실습폴더에 **23) comp_table_L6.qza** 폴더가 새로 생성된 것을 확인하세요.

```sh
qiime composition ancom \
--i-table comp_table_L6.qza \
--m-metadata-file sample-metadata_ata.tsv \
--m-metadata-column Vegetation \
--o-visualization ancom-vegetation_L6.qzv
```
- 현재 실습폴더에 **24) ancom-vegetation_L2.qzv** 폴더가 새로 생성된 것을 확인하세요.
- 해당 .qzv 파일을 [QIIME2view](https://view.qiime2.org) 에서 열어보세요.   

</br></br>

# 9. Exporting 

QIIME2 이외에 다른 툴을 이용하여 추가 분석하고 싶을 때, qza 포멧이 아닌 biom 또는 txt 포멧이 필요할 수 있습니다.

> ### *.qza 파일(count data)을 *.biom 포멧으로 Export 하기
    - 필요한 파일
        - FeatureTable[Frequency] : 09) filtered_100_table.qza -> features 의 read count 
        - FeatureData[Taxonomy] : 16) silva138_99_taxonomy.qza


- Feature table export
```
qiime tools export \
--input-path filtered_100_table.qza \
--output-path exported_biom
```

- Taxonomy export

```
qiime tools export \
--input-path silva138_99_taxonomy.qza \
--output-path exported_biom
```
- 현재폴더에 25) exported_biom 폴더 생성된 것을 확인 후, exported_biom/taxonomy.tsv 를 ./biom-taxonomy.tsv 으로 복사본 만들기
```
cp ./exported_biom/taxonomy.tsv ./biom-taxonomy.tsv
```

- 현재폴더(./)에 copy 된 26) biom-taxonomy.tsv 파일을 엑셀이나 메모장으로 열기

Feature ID | Taxon | Confidence  
---|---|---

로 되어 있는 Header 를

#OTUID | taxonomy | confidence  
---|---|---

로 변경한 후 biom-taxonomy.tsv 로 그대로 저장(덮어쓰기)하고 아래 명령어 실행


```
biom add-metadata \
-i ./exported/feature-table.biom \
-o ./table-with-taxonomy.biom \
--observation-metadata-fp ./biom-taxonomy.tsv \
--sc-separated taxonomy
```

- filtered_100_table.qza 이 27) table-with-taxonomy.biom 로 변경 완료.
- QIIME2 에서는 feature table 과 taxonomy 파일이 별도로 존재했지만, biom 파일은 이 두가지 파일이 biom 파일 하나에 모두 들어가 있습니다.
</br></br>


> ### *.qza 파일(read count)을 relative abundance.qza 파일로 변경한 후, *.tsv 포멧으로 Export 하기

[MaAsLin](https://huttenhower.sph.harvard.edu/maaslin/) 이나 [LEfSe](https://github.com/biobakery/biobakery/wiki/lefse#3-visualization) 등의 다른 프로그램으로 분석하고자 할때, read count 데이터가 아닌 relative abundance (proportional data) 데이터가 사용됩니다. count 데이터를 relative abundance 로 변경 후 export 진행합니다. </br></br> 
또한 qza 파일이나 biom 파일 등은 데이터를 바로 열어 확인하기 어려우나, tsv 파일포멧으로 변환하면 엑셀이나 텍스트 프로그램에서 데이터를 눈으로 바로 확인 가능합니다.

- *.qza -> *.biom -> *.tsv 단계로 진행함.
- 위에서 biom 파일 exporting 할 때와 같이 ASV 파일을 export 할 수도 있지만, taxonomy level 별로 별도의 파일로도 만들 수 있음

- 원하는 taxa level (예, Phylum: L2)로 collapse
```
qiime taxa collapse \
--i-table filtered_100_table.qza \
--i-taxonomy silva138_99_taxonomy.qza \
--p-level 2 \
--o-collapsed-table L2-table.qza
```
- output 파일 28) L2-table.qza 생성

- read count 데이터인 Featuretable[Frequency] 를 relative frequency 데이터인 FeatureTable[RelativeFrequency] 로 변경
```
qiime feature-table relative-frequency \
--i-table L2-table.qza \
--o-relative-frequency-table R-L2-table.qza
```
- 생성된 output 파일 29) R-L2-table.qza 를 biom 파일로 export
- taxa level 로 collapse 된 파일은 taxonomy 파일을 별도로 export 할 필요 없음

```
qiime tools export \
--input-path R-L2-table.qza \
--output-path R_L2
```
- 생성된 R_L2 폴더 내에 30) feature-table.biom 생성 확인
```
cd ./R_L2
```
- feature-table.biom 을 .tsv 파일로 변경
```
biom convert -i feature-table.biom -o R-L2-table.tsv --to-tsv
```
- output 파일 31) R-L2-table.tsv 파일을 엑셀이나 메모장으로 열어서 확인합니다. 
- L2 뿐만아니라 L3 (class), L4 (order), L5 (family), L6 (genus), L7 (species) 까지 같은 방식으로 Exporting 가능합니다.
- 이 taxa level 별 relative abundance 파일은 그룹별 abundance 비교 box plot 을 새로 그리는 등 다양하게 활용 가능합니다. 


</br></br>

## 수고하셨습니다! </br>
### 본 스크립트는 잘 보관하셔서 차후 QIIME2 분석하실 때 요긴하게 사용하시길 바랍니다~:poop::sweat_smile::grin::heart:

---
본 메뉴얼의 저작권은 [QIIME2](https://docs.qiime2.org/2019.4/citation/) 에 있습니다. <br>
강의자료 제작: <span style="color:blue">강북삼성병원 **김한나**</span> (문의: 147942@hanmail.net) 
---

