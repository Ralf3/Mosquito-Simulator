# Run a Simulation
## Change the control.jl

Please open the control file control.jl and have a look at it. Control.jl controls the simulation and contains all configuration options. Before we start a simulation, you can make a small adjustment. We want to add 20 infected birds to the simulation region at step 20, i.e. 20 days after the start of the simulation.

- infstepIB=20
- addIB=20

Then save control.jl

### Start the Simulation
**include("control.jl")**

Normally the simulation starts very fast. However, the first time some Julia packages have to be compiled, which may take a few minutes. Please be patient. You can see that the simulation is running by the output of the simulation in form of a table. This is in my eyes better than a progress bar.

Wenn die Simulation beendet wurde wird ein Bild **show_resb()** erscheinen. Es zeigt die Dynamik des WestNil Fieber bei den Vögeln. Deutlich ist der Pik um den Step 20 zu erkennen. Dort wurden 20 infizierte Vögel eingeführt.

<p><img src="doc/picts/show_resb.png",title="Bird dynamics"/></p>