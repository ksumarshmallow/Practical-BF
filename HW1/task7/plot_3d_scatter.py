import numpy as np
import plotly.graph_objects as go

def plot_3d_scatter(x: np.ndarray, y: np.ndarray, z: np.ndarray,
                    x_title: str, y_title: str, z_title: str,
                    c: np.ndarray, plot_name: str):
    
    fig = go.Figure(data=[go.Scatter3d(
        x=x, y=y, z=z,
        mode='markers',
        marker=dict(
            size=5,
            color=c,          
            colorscale='Viridis',   
            colorbar=dict(title='POS (chr15, hg38)')
        )
    )])

    # Обновление осей
    fig.update_layout(
        scene=dict(
            xaxis_title=x_title,
            yaxis_title=y_title,
            zaxis_title=z_title
        ),
        title=plot_name
    )

    # Меняем размера графика
    fig.update_layout(
        width=1000, 
        height=800
    )

    return fig