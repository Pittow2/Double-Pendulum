角度 $\alpha,\beta$

角速度 $\alpha_0,\alpha_1$

质量 $m_0,m_1$

长度 $r_0,r_1$

重力加速度 $m_0:(0,g_0),m_1:(g_x,g_y)$

杆 $1$ 拉力 $x$

$m_0$ 加速度（向心，加速）： $\alpha_0^2 r_0,-g_0\sin\alpha-\frac{x}{m_0}\sin(\alpha-\beta)$

$m_1$ 加速度（向心，加速）： $\frac{x}{m_1}-g_x\sin\beta-g_y\cos\beta,g_x\cos\beta-g_y\sin\beta$

杆 $1$ 向心加速度： 

$$
\frac{x}{m_1}-g_x\sin\beta-g_y\cos\beta-\alpha_0^2 r_0\cos(\alpha-\beta)+\sin(\alpha-\beta)(g_0\sin\alpha+\frac{x}{m_0}\sin(\alpha-\beta))=\\
x(\frac{1}{m_1}+\frac{\sin^2(\alpha-\beta)}{m_0})-g_x\sin\beta-g_y\cos\beta-\alpha_0^2 r_0\cos(\alpha-\beta)+g_0\sin\alpha\sin(\alpha-\beta)=
\alpha_1^2 r_1
$$

解得：

$$
x=\frac{g_x\sin\beta+g_y\cos\beta+\alpha_0^2 r_0\cos(\alpha-\beta)+\alpha_1^2 r_1-g_0\sin\alpha\sin(\alpha-\beta)}{\frac{1}{m_1}+\frac{\sin^2(\alpha-\beta)}{m_0}}
$$

角加速度：
$$
a_0=-\frac{g_0\sin\alpha+\frac{x}{m_0}\sin(\alpha-\beta)}{r_0}
\\
a_1=\frac{g_x\cos\beta-g_y\sin\beta+\alpha_0^2 r_0\sin(\alpha-\beta)-a_0 r_0\cos(\alpha-\beta)}{r_1}
$$
一些简化：
$$
\frac{x}{m_0}=\frac{g_x\sin\beta+g_y\cos\beta+\alpha_0^2 r_0\cos(\alpha-\beta)+\alpha_1^2 r_1-g_0\sin\alpha\sin(\alpha-\beta)}{\frac{m_0}{m_1}+\sin^2(\alpha-\beta)}
$$
若 $c=\alpha_0,d=\alpha_1,e=\frac{m_0}{m_1},f=\frac{r_1}{r_0},s_0=\sin\alpha,s_1=\sin\beta,c_1=\cos\beta,s_t=\sin(\alpha-\beta),c_t=\cos(\alpha-\beta)$

且 $r_0=1$
$$
A_0=g_0 s_0+\frac{s_t(rc^2 c_t+rfd^2-g_0 s_0 s_t+g_x s_1+g_y c_1)}{e+s_t^2}=g_0 s_0+\frac{s_t(c^2 c_t+fd^2-g_0 s_0 s_t+g_x s_1+g_y c_1)}{e+s_t^2}
\\
a_0=-\frac{A_0}{r}=-A_0
\\
a_1=\frac{g_x c_1-g_y s_1+rc^2 s_t+A_0 c_t}{rf}=\frac{g_x c_1-g_y s_1+c^2 s_t+A_0 c_t}{f}
$$


