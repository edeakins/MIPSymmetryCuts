\Problem name: ClpDefaultName

Minimize
obj: x0 + x1 + x2 + x3 + x4 + x5
Subject To
cons0:  x0 + x1 + x2 + x3 >= 1
cons1:  x0 + x1 + x2 + x4 >= 1
cons2:  x0 + x1 + x2 + x5 >= 1
cons3:  2 x0 + 2 x3 + x4 + x5 >= 2
cons4:  2 x1 + x3 + 2 x4 + x5 >= 2
cons5:  2 x2 + x3 + x4 + 2 x5 >= 2
cons6:  -2 x0 -3 x1 -3 x2 -3 x3 -3 x4 -3 x5 <= -6
cons7:  -4 x0 -4 x1 -4 x2 -4 x3 -4 x4 -4 x5 <= -8
cons8:  -8 x0 -6 x1 -12 x2 -12 x3 -12 x4 -15 x5 <= -24
cons9:  -2 x0 -2 x2 -2 x3 -2 x4 -3 x5 <= -4
Bounds
 0 <= x0 <= 1
 0 <= x1 <= 1
 0 <= x2 <= 1
 0 <= x3 <= 2
 0 <= x4 <= 2
 0 <= x5 <= 2
Integers
x0 x1 x2 x3 x4 x5 
End
