import numpy as np
jj = np.empty(16)
ii = np.empty(16)
m = 0
for j in [0,1]:
  for i in [0,1]:
    ii[m] = i
    jj[m] = j
    m = m + 1
for j in [-1,0]:
  for i in [0,1]:
    ii[m] = i
    jj[m] = j
    m = m + 1
for j in [0,1]:
  for i in [-1,0]:
    ii[m] = i
    jj[m] = j
    m = m + 1
for j in [-1,0]:
  for i in [-1,0]:
    ii[m] = i
    jj[m] = j
    m = m + 1

num = np.empty((3, 3))
for j in range(-1, 2):
  for i in range(-1, 2):
    num[i+1, j+1] = 0
    for m in range(16):
      if ii[m] == i and jj[m] == j:
        num[i+1,j+1] += 1


for j in range(-1, 2):
  for i in range(-1, 2):
    print(f"({i}, {j}): {num[i+1, j+1]}")


# print("----")
# for j in range(-1, 2):
#   for i in range(-1, 2):
#     if num[i+1, j+1]==8.0:
#       print(f"({i}, {j}, {k})")