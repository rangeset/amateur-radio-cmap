> 中文编写 英语不好 有英语好的可以翻译下提个PR 感谢

其实原理很简单

首先我们有一个球面 $ S_0 $ 半径为 $ R $

然后球上一点 $ P_0(\theta _0, \varphi _0) $ 当作初始位置

然后 我们从这个点作另一个球面 $ S $ 两球面交圆为 $ C $

该圆上的点到 $ S_0 $ 距离都为 $ S $的半径 $ r $

容易得到我们的实际绘图的距离中点的距离为 $ \rho = arcsin(r) $

现在我们将过点$ P_0 $ 与 $ S_0 $ 相切 且z轴正半轴相交的线作为我们绘图的y轴

设 $ \theta_p, \varphi_p $ 分别作为 $ C $ 上点的纬度和经度

容易得到

$$
r = \sqrt{sin^2(\theta_p - \theta_0) + sin^2(\varphi_p - \varphi_0)} \\
\rho = arcsin(\sqrt{sin^2(\theta_p - \theta_0) + sin^2(\varphi_p - \varphi_0)}) \\
\theta = arctan(\frac{sin(\varphi_p - \varphi_0)}{sin(\theta_p - \theta_0)})
$$

代入求解

$$
\begin{aligned}
x &= \rho \cdot cos(\theta) \\
&= arcsin(\sqrt{sin^2(\theta_p - \theta_0) + sin^2(\varphi_p - \varphi_0)}) \cdot cos(arctan(\frac{sin(\varphi_p - \varphi_0)}{sin(\theta_p - \theta_0)})) \\
&= \frac{arcsin(\sqrt{sin^2(\theta_p - \theta_0) + sin^2(\varphi_p - \varphi_0)})}{\sqrt{1+(\frac{sin(\varphi_p - \varphi_0)}{sin(\theta_p - \theta_0)})^2}} \\
&= \frac{arcsin(\sqrt{sin^2(\theta_p - \theta_0) + sin^2(\varphi_p - \varphi_0)})}{\sqrt{sin(\varphi_p - \varphi_0)^2 + sin(\theta_p - \theta_0)^2}} \cdot sin(\theta_p - \theta_0) \\
&= \rho \cdot sin(\theta - \theta_0) \\
&= \frac{arcsin(\sqrt{sin^2(\theta_p - \theta_0) + sin^2(\varphi_p - \varphi_0)})}{\sqrt{sin(\varphi_p - \varphi_0)^2 + sin(\theta_p - \theta_0)^2}} \cdot sin(\varphi_p - \varphi_0) \\
\end{aligned}
$$

得到最后的方程

$$
\begin{cases}
x = \frac{arcsin(\sqrt{sin^2(\theta_p - \theta_0) + sin^2(\varphi_p - \varphi_0)})}{\sqrt{sin(\varphi_p - \varphi_0)^2 + sin(\theta_p - \theta_0)^2}} \cdot sin(\theta_p - \theta_0) \\
y = \frac{arcsin(\sqrt{sin^2(\theta_p - \theta_0) + sin^2(\varphi_p - \varphi_0)})}{\sqrt{sin(\varphi_p - \varphi_0)^2 + sin(\theta_p - \theta_0)^2}} \cdot sin(\varphi_p - \varphi_0)
\end{cases}
$$