
import matplotlib.pyplot as plt
mald = True

thrust_values = []
cc_pressure_values = []
timestamps = []

plt.ion()
fig1 = plt.figure(animated=True)

# fig1, axs = plt.subplots(1,2)

ax1 = fig1.add_subplot(1,2,1)
ax2 = fig1.add_subplot(1,2,2)

# ax1 = axs[0]
# ax2 = axs[1]

plot, = ax1.plot(timestamps, thrust_values)
plot2, = ax2.plot(timestamps, cc_pressure_values)



ax1.set_title('Thrust')
ax1.set_xlim(0,100)

ax2.set_title('CC Pressure')
ax2.set_xlim(0,100)

fig1.canvas.draw()


plt.show(block=False)

time = 1

while mald:
    time += 1
    timestamps.append(time)
    thrust_values.append(1)
    cc_pressure_values.append(2)

    plot.set_data(timestamps, thrust_values)
    plot2.set_data(timestamps, cc_pressure_values)

    # ax1.redraw_in_frame()
    # ax2.redraw_in_frame()

    fig1.canvas.draw()
    fig1.canvas.flush_events()  

    print('Mald ' + str(time))


