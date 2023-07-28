# MBET

## Build & Run

```shell
mkdir build && cd build
cmake ..
make MBET MBETM

./MBET ../data/youtube-groupmemberships.txt 1
./MBETM ../data/github.txt 1
```

To output results of the toy graph in the paper:

```shell
make MBETR
./MBETR ../data/test.txt 1
```
