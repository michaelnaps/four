## Fourier Series and Animations

**UNDER DEVELOPMENT**

This small Python library was inspired by the **3Brown1Blue** video linked [here](https://www.youtube.com/watch?v=-qgreAUpPwM). After watching this video I took some time to remind myself of the structure of a Fourier series (which I, at one point, taught in ungrad.) before building the classes in the package **FOUR**.

The idea of a 2-D Fourier series is relatively simple. The general theory states that any function, $f$, which maps between two domains such that $f : \mathbb{X} \rightarrow \mathbb{Y}$ can be represented by an infinite sum of $\cos(x)$ and $\sin(x)$ terms.

$$
    f(x) = \sum_{k=0}^\infty \left[ A_k \sin(\frac{\pi k x}{L}) + B_k \cos(\frac{\pi k x}{L}) \right]
$$

Where $x \in \mathbb{X}$ and $f(x) \in \mathbb{Y}$.


___

### Imprinting Temporal-dynamics on to a Grid

In this quick example I have imprinted the motion of the duffing oscillator (as begun from a single set point) on to a grid. Each space can be made *visible* of *invisible* dependent on whether the point is present in that cell. That is,

$$
    c_{i,j} = \begin{cases}
        1 & \text{if $x(t)$ is present} \\
        0 & \text{otherwise}
    \end{cases}
$$

Where $c_{i,j}$ is the binary value for the $j$-th cell in the $i$-th row. The value of each cell was saved for every point of the simulation to be used in the calculation of the Fourier transform.

Then, a Fourier transform was performed for each cell in the grid, resulting in the figure below...

<p align="center">
    <img src=duffing/grid.gif width=300 />
</p>


___

### Abby and Michael Simulated Trace

**Note:** I did this purely as a side project. Animations using Fourier transforms are not particularly difficult, but the exercise proved a fun and interesting re-introduction to a topic I am quite enamored with presently. For the purposes of the simulation I chose $N=50$ in my least-squares (LS) regression. The discrete fourier transform (DFT) is also a method in this library, but was not used here as I wanted to filter out some of the hard-edges from my trace.

For Abby's birthday I traced an image using the tools in the *drawdata* folder of this repository to form a reference data set, and used the *FOUR*-package to create an animation of the lines drawn. The first image shows the reference...

<p align="center">
    <img src=.archive/draw/images/abby_michael.png width=300 />
</p>

And the second shows the resulting animation (plus I gave myself sunglasses)...

<p align="center">
    <img src=abby_pkg/recording.gif width=450 />
</p>

NOTE: Simulation is very smooth when not formatted as a *.gif* file.
