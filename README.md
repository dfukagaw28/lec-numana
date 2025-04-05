# 同志社大学文化情報学部『数値解析』講義資料（2025年度春学期）

## How to build

### Build the book with podman

Build a container image:

```
(in this directory)

$ podman build -t myjb .
```

Create and run a docker container:

```
$ podman run --rm -it -v $PWD:/work -w /work myjb jupyter-book build docs --all
```
