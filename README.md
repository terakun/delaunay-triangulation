# Delaunay三角形分割

`input.dat`:
```
0 0 0
1 1 0
2 0.5 0.5
3 0 1
4 1 1
```
![points](https://i.imgur.com/p0ztzWe.png)

```
$ ./main input.dat
```

`plot.txt`:
```
0 0
1 0
0.5 0.5
0 0

0.5 0.5
0 1
0 0
0.5 0.5

1 1
0.5 0.5
1 0
1 1

1 1
0 1
0.5 0.5
1 1
```
![triangle](https://i.imgur.com/vmq83Qf.png)

`graph.txt`:
```
5
3 1 2 3
3 0 2 4
4 0 1 3 4
3 0 2 4
3 1 2 3
```
