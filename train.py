

import matplotlib.pyplot as plt

def sfd_bmd(offset):
  P = 400; # Total weight of train [N]
  L = 1260; # Length of bridge
  x = range(0, L + 1); # x-axis
  x_train = array_add([52, 228, 392, 568, 732, 908], offset)
  P_train = array_mul([1, 1, 1, 1, 1, 1],P/6)

  F_Ay, F_By = None, None

  M_A = 0
  for loadingNumber in range(len(x_train)):
    M_A += x_train[loadingNumber] * P_train[loadingNumber]

  # All the loadings are -distance*weight but when you solve the equation (move all loadings to other side) they become positive. 
  F_By = M_A / 1200
  F_Ay = P - F_By


  loadings = [F_Ay, F_By] + array_mul(P_train, -1) 
  loading_positions = [0, L] + x_train
  SFD = []
  for pos in x:
    SFD.append(add_loadings(pos, loadings, loading_positions))
  # f = open('sfd_output.txt', 'w')
  # for value in SFD:
  #   # print(f'writing{value}')
  #   f.write(f"{value}\n")
  # f.close()
  

  BMD = [0] * len(SFD)
  r_start, r_end = 0, 0
  while r_end < len(SFD):
      # print(r_end)
      if SFD[r_start] == SFD[r_end]:
        r_end += 1
      else:
        # At this point the start will be the value that it was previously but the end will be the value of the next rectangle
        change_in_bending = (r_end - r_start) * SFD[r_start]
        add_bending_change(BMD, change_in_bending, r_start, r_end)
        propogate_moment(BMD, r_end, BMD[r_end])
        if (r_end < len(SFD)):
          r_start = r_end # Start the next rectangle at that point
      #print(i)

  # add final change

  # print("DEBUG: last change", r_start,r_end)
  # change_in_bending = (r_end - r_start) * SFD[r_start]

  # add_bending_change(BMD, change_in_bending, r_start, r_end)
  # propogate_moment(BMD, r_end, BMD[r_end])
  # # print(SFD, BMD)
  # with open('bmd_output.txt', 'w') as f:
  #   for value in BMD:
  #     f.write(f"{value}\n")
  # f.close()

  return max(SFD), max(BMD)

def add_loadings(cur_pos, loadings, loading_positions):
  loading = 0
  for i in range(len(loadings)):
    loading_pos = loading_positions[i]
    if (cur_pos >= loading_pos):
      loading += loadings[i]
  return loading

def array_mul(arr, constant):
  for i in range(len(arr)):
    arr[i] *= constant
  return arr
def array_add(arr, constant):
  for i in range(len(arr)):
    arr[i] += constant
  return arr

def add_bending_change(BMD, change_in_bending, r_start, r_end):
  sum = 0 
  slope =  change_in_bending / (r_end - r_start)
  for i in range(r_start, r_end + 1):
    BMD[i] += slope * (i - r_start)

def propogate_moment(BMD, start, val):
    '''Propogate the current moment value across the rest of the positions, this will make it so that you're adding some dy to the existing moment and not zero'''
    for i in range(start, len(BMD)):
      BMD[i] = val
    return BMD


def sfd_bmd_env():
  sfds = []
  bmds =[]
  for offset in range(-52, 292 + 1):
    sfd_max, bmd_max = sfd_bmd(offset)
    sfds.append(sfd_max)
    bmds.append(bmd_max)

  
  plt.figure(figsize=(12, 6))

  plt.subplot(1, 2, 1)
  plt.plot(range(0, 292 + 1), sfds, label='SFD Max')
  plt.xlabel('Offset')
  plt.ylabel('SFD Max')
  plt.title('SFD Max vs Offset')
  plt.legend()

  plt.subplot(1, 2, 2)
  plt.plot(range(0, 292 + 1), bmds, label='BMD Max')
  plt.xlabel('Offset')
  plt.ylabel('BMD Max')
  plt.title('BMD Max vs Offset')
  plt.legend()

  plt.tight_layout()
  plt.show()
sfd_bmd_env()


