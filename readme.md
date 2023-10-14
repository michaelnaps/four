## Fourier Series and Animations

**UNDER DEVELOPMENT**

This small Python library was inspired by the **3Brown1Blue** video linked [here](https://www.youtube.com/watch?v=-qgreAUpPwM). After watching this video I took some time to remind myself of the structure of a Fourier series (which I, at one point, taught in ungrad.) before building the classes in the package **FOUR**.

The idea of a 2-D Fourier series is relatively simple. The general theory states that any function, $f$, which maps between two domains such that $f : \mathbb{X} \rightarrow \mathbb{Y}$ can be represented by an infinite sum of $\cos(x)$ and $\sin(x)$ terms.

$$
    f(x) = \sum_{k=0}^\infty \left[ A_k \sin(\frac{\pi k x}{L}) + B_k \cos(\frac{\pi k x}{L}) \right]
$$

Where $x \in \mathbb{X}$ and $f(x) \in \mathbb{Y}$.

    References:
        A. Gilat and V. Subramaniam, Numerical methods for engineers and
            scientists: an introduction with applications using matlab, 3. ed.
            Hoboken, NJ: Wiley, 2014.

___

### Imprinting Temporal-dynamics onto a Grid

In this quick example I have imprinted the motion of the duffing oscillator (as begun from a single set point) onto a grid. Each space can be made *visible* of *invisible* dependent on whether the point is present in that cell. That is,

$$
    c_{i,j} = \begin{cases}
        1 & \text{if $x(t)$ is present} \\
        0 & \text{otherwise}
    \end{cases}
$$

Where $c_{i,j}$ is the binary value for the $j$-th cell in the $i$-th row. The binary grid was saved at each step of the simulation and used to calculate the Fourier transform. This resulted in a playback that looks like...

<p align="center">
    <img src=duffing/grid.gif width=300 />
</p>


___

### Abby and Michael Simulated Trace

This started as a fun side project while working at iRobot in the Summer of 2023. Animations using Fourier transforms are not particularly difficult, but the exercise proved an interesting re-introduction to a topic I was enamored with at the time. 

For Abby's birthday I traced an image using the tools in the *drawdata* folder of this repository to form a reference data set, and used the *FOUR*-package to create an animation of the lines drawn. For the purposes of the simulation I chose $N=50$ in my least-squares (LS) regression. The discrete fourier transform (DFT) is also a method in this library, but was not used here as I wanted to filter out some of the hard-edges from my trace. The first image shows the reference...

<p align="center">
    <img src=.archive/draw/images/abby_michael.png width=300 />
</p>

And the second shows the resulting animation (plus I gave myself sunglasses)...

<p align="center">
    <img src=abby_pkg/recording.gif width=450 />
</p>

NOTE: Simulation is very smooth when not formatted as a *.gif* file.
