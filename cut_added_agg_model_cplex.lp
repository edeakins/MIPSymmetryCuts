\ENCODING=ISO-8859-1
\Problem name: 

Minimize
 obj1: x0 + x1 + x2 + x6 + x7 + x8
Subject To
 R0#0:   x0 + x1 + x2 + x6 >= 1
 R1#1:   x0 + x1 + x2 + x7 >= 1
 R2#2:   x0 + x1 + x2 + x8 >= 1
 R6#3:   2 x0 + 2 x6 + x7 + x8 >= 2
 R7#4:   2 x1 + x6 + 2 x7 + x8 >= 2
 R8#5:   2 x2 + x6 + x7 + 2 x8 >= 2
Bounds
 0 <= x0 <= 1
 0 <= x1 <= 1
 0 <= x2 <= 1
 0 <= x6 <= 2
 0 <= x7 <= 2
 0 <= x8 <= 2
End
