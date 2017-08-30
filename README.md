# naive-binary-gemm

The results of `cblas_sgemm` are two baselines:

1. Performance(speed) benchmark: don't care result precision of `cblas_sgemm`;  
2. Precision benckmark: validate the result of naive binary GEMM.

## Build

```shell
./make.sh
```

## Usage

```shell
./naive-binary-gemm
```
