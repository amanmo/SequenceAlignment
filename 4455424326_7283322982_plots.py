import matplotlib.pyplot as plt

def plot_spacetime():

    index = [128,256,512,768,1024,1536,2048,3072,4096,6144]

    time_taken_DP = [0.012668848037719727,0.04834389686584473,0.18639302253723145,0.3816568851470947,0.7581849098205566,1.516230821609497,3.041546106338501,6.148684024810791,12.171642780303955,24.29126000404358]
    time_taken_DC = [0.02377486228942871,0.09153032302856445,0.3621518611907959,0.7546999454498291,1.4760360717773438,2.9603331089019775,6.148335933685303,12.037397146224976,24.1245219707489,49.14292097091675]

    space_taken_DP = [352,864,2400,3728,4992,8976,13920,26304,52608,106176]
    space_taken_DC = [32,80,128,190,256,376,528, 848, 966,1126]

    plt.plot(index, time_taken_DP, label='Dynamic Programming')
    plt.plot(index, time_taken_DC, label='Divide and Conquer')
    plt.xlabel('Input Size')
    plt.ylabel('Time Taken (seconds)')
    plt.title('Time Complexity')
    plt.legend()
    plt.savefig('Time Complexity.png')
    plt.close()

    plt.plot(index, space_taken_DP, label='Dynamic Programming')
    plt.plot(index, space_taken_DC, label='Divide and Conquer')
    plt.xlabel('Input Size')
    plt.ylabel('Space Consumed (kilobytes)')
    plt.title('Space Complexity')
    plt.legend()
    plt.savefig('Space Complexity.png')
    plt.close()

plot_spacetime()