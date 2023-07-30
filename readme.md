## Fourier Series and Animations

**UNDER DEVELOPMENT**

This small Python library was inspired by the **3Brown1Blue** video linked [here](https://www.youtube.com/watch?v=-qgreAUpPwM). After watching this video I took some time to remind myself of the structure of a Fourier series (which I, at one point, taught in ungrad.) before building the classes in the package **FOUR**.

The idea of a 2-D Fourier series is relatively simple. The general theory states that any function, $f$, which maps between two domains such that $f : \mathbb{X} \rightarrow \mathbb{Y}$, then the function can be represented by an infinite sum of $\cos(x)$ and $\sin(x)$ terms.

$$
    f(x) = \sum_{k=0}^\infty \left[ A_k \sin(\frac{\pi k x}{L}) + B_k \cos(\frac{\pi k x}{L}) \right]
$$

Where $x \in \mathbb{X}$ and $f(x) \in \mathbb{Y}$.


___

### Abby and Michael Simulated Trace

<video width="320" height="240" controls>
  <source src="abby_pkg/recording.mp4" type="video/mp4">
</video>