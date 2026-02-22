Overlap test - two spheres sharing a region
c
c Cell cards
c
1  1  -2.7   -1          $ sphere at origin r=5, material 1
2  2  -8.0   -2          $ sphere at (3,0,0) r=5, material 2 - overlaps cell 1
3  3  -1.0   -3  1  2    $ box shell: inside box but outside both spheres
99 0          3           $ graveyard (outside box)

c Surface cards
c
1  s   0  0  0  5        $ sphere at origin, radius 5
2  s   0  0  0  3        $ sphere at (3,0,0), radius 5
3  rpp -10 15 -10 10 -10 10  $ bounding box

c Data cards
m1  26056.80c 1.0        $ iron
m2  74184.80c 1.0        $ tungsten
m3  1001.80c  1.0        $ hydrogen
