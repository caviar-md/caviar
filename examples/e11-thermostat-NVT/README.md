


This example contains two stages thermostat. First it uses Berendesen to reach the correct level of Temperature. Then it removes Berendesen and adds Nose-Hoover to create a correct canonical fluctuation for the temperature. Using only Nose-Hoover may make problem.
The exported file 'o_temperature.txt' contains the temperature and it can be plotted via gnuplot command:
  plot 'o_temperature.txt' u 1:3 w l


It is developed according to the following reference:
" Berendsen and Nose-Hoover thermostats "
  Victor Rühle "
" August 8, 2007 "
