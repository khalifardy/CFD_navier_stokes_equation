import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import seaborn as sns
sns.set()

N_POINTS = 41
DOMAIN_SIZE = 1.0
N_ITERATIONS = 500
TIME_STEP_LENGTH = 0.001
KINEMATIC_VISCOSITY = 0.1
HORIZONTSAL_VELOCITY_TOP = 1.0
DENSITY = 1
STABILITY_SAFETY_FACTOR = 0.5

N_PRESSURE_POISSON_ITERATIONS = 50
def central_difference_x(f,element_length):
    diff = np.zeros_like(f)
    diff[1:-1, 1:-1] = (
        f[1:-1, 2:]
        -
        f[1:-1, 0:-2]
    )/(
        2*element_length
    )
    return diff

def central_difference_y(f,element_length):
    diff = np.zeros_like(f)
    diff[1:-1, 1:-1] = (
        f[2: , 1:-1]
        -
        f[0:-2, 1:-1]
    )/(
        2*element_length
    )
    return diff

def laplace(f,element_length):
    diff = np.zeros_like(f)
    diff[1:-1, 1:-1] = (
        f[1:-1 , 0:-2]
        +
        f[0:-2, 1:-1]
        -
        4
        *
        f[1:-1, 1:-1]
        +
        f[1:-1, 2:]
        +
        f[2:,1:-1]
    )/(
        element_length**2
    )
    return diff


def main():
    element_length = DOMAIN_SIZE / (N_POINTS-1)
    x = np.linspace(0.0,DOMAIN_SIZE,N_POINTS)
    y = np.linspace(0.0, DOMAIN_SIZE, N_POINTS)

    X,Y = np.meshgrid(x,y)
    u_prev = np.zeros_like(X)
    v_prev = np.zeros_like(X)
    p_prev = np.zeros_like(X)

    maximum_possible_time_step_length = (
        0.5 * element_length**2/KINEMATIC_VISCOSITY
    )

    if TIME_STEP_LENGTH > STABILITY_SAFETY_FACTOR * maximum_possible_time_step_length:
        raise RuntimeError('Stabilitas tidak di garansi')

    for _ in tqdm(range(N_ITERATIONS)):
        d_u_prev__d_x = central_difference_x(u_prev,element_length)
        d_u_prev__d_y = central_difference_y(u_prev,element_length)
        d_v_prev__d_x = central_difference_x(v_prev,element_length)
        d_v_prev__d_y = central_difference_y(v_prev,element_length)
        laplace__u_prev = laplace(u_prev,element_length)
        laplace__v_prev = laplace(v_prev,element_length)

        u_tent=(
            u_prev
            +
            TIME_STEP_LENGTH *(
                -
                (
                    u_prev * d_u_prev__d_x
                    +
                    v_prev*d_u_prev__d_y
                )
                +
                KINEMATIC_VISCOSITY *laplace__u_prev
            
            )
        
        )

        v_tent=(
            v_prev
            +
            TIME_STEP_LENGTH *(
                -
                (
                    u_prev * d_v_prev__d_x
                    +
                    v_prev*d_v_prev__d_y
                )
                +
                KINEMATIC_VISCOSITY *laplace__v_prev
            
            )
        
        )

        u_tent[0, :] = 0
        u_tent[:, 0] = 0
        u_tent[:, -1] = 0
        u_tent[-1, :]=HORIZONTSAL_VELOCITY_TOP

        v_tent[0, :] = 0
        v_tent[:, 0] = 0
        v_tent[:, -1] = 0
        v_tent[-1, :]= 0

        d_u_tent__d_x = central_difference_x(u_tent,element_length)
        d_v_tent__d_y = central_difference_y(v_tent,element_length)

        rhs = (
            DENSITY / TIME_STEP_LENGTH
            *
            (
                d_u_tent__d_x
                +
                d_v_tent__d_y
            )
        )

        for _ in range(N_PRESSURE_POISSON_ITERATIONS):
            p_next = np.zeros_like(p_prev)
            p_next[1:-1,1:-1] = 1/4 *(
                +
                p_prev[1:-1,0:-2]
                +
                p_prev[0:-2, 1:-1]
                +
                p_prev[1:-1, 2: ]
                +
                p_prev[2:, 1:-1]
                -
                element_length**2
                *rhs[1:-1,1:-1]
            )
        

            p_next[:, -1] = p_next[:, -2]
            p_next[0, :] = p_next[1, :]
            p_next[:, 0] = p_next[:, 1]
            p_next[-1, :] = 0
            
            p_prev = p_next
        
        d_p_next__d_x = central_difference_x(p_next,element_length)
        d_p_next__d_y = central_difference_y(p_next,element_length)

        u_next = (
            u_tent
            -
            TIME_STEP_LENGTH/DENSITY
            *
            d_p_next__d_x
        )

        v_next = (
            v_tent
            -
            TIME_STEP_LENGTH/DENSITY
            *
            d_p_next__d_y
        )

        u_next[0, :] = 0
        u_next[:, 0] = 0
        u_next[:, -1] = 0
        u_next[-1, :]=HORIZONTSAL_VELOCITY_TOP

        v_next[0, :] = 0
        v_next[:, 0] = 0
        v_next[:, -1] = 0
        v_next[-1, :]= 0

        u_prev = u_next
        v_prev = v_next
        p_prev = p_next
    
    plt.figure()
    plt.contourf(X, Y, p_next)
    plt.colorbar()

    plt.quiver(X,Y, u_next,v_next, color='black')
    plt.streamplot(X,Y, u_next, v_next, color='black')
    plt.show()



if __name__ =='__main__':
    main()