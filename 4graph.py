import matplotlib.pyplot as plt

# Data from your results
positions = list(range(131))  # From 0 mm to 130 mm

# Analytical solution
analytical_temps = [
    600.000, 599.540, 599.080, 598.620, 598.160, 597.700, 597.240, 596.780, 596.320, 595.860,
    595.400, 593.586, 591.771, 589.955, 588.139, 586.323, 584.507, 582.692, 580.876, 579.060,
    577.244, 575.428, 573.613, 571.797, 569.981, 568.165, 566.349, 564.534, 562.718, 560.902,
    559.086, 557.271, 555.455, 553.639, 551.823, 550.007, 548.192, 546.376, 544.560, 542.744,
    540.928, 539.113, 537.297, 535.481, 533.665, 531.849, 530.034, 528.218, 526.402, 524.586,
    522.771, 520.955, 519.139, 517.323, 515.507, 513.692, 511.876, 510.060, 508.244, 506.428,
    504.613, 502.797, 500.981, 499.165, 497.349, 495.534, 493.718, 491.902, 490.086, 488.271,
    486.455, 484.639, 482.823, 481.007, 479.192, 477.376, 475.560, 473.744, 471.928, 470.113,
    468.297, 466.481, 464.665, 462.849, 461.034, 459.218, 457.402, 455.586, 453.771, 451.955,
    450.139, 448.323, 446.507, 444.692, 442.876, 441.060, 439.244, 437.428, 435.613, 433.797,
    431.981, 430.165, 428.349, 426.534, 424.718, 422.902, 421.086, 419.271, 417.455, 415.639,
    413.823, 413.110, 412.420, 411.730, 411.040, 410.350, 409.660, 408.970, 408.280, 407.590,
    406.900, 406.210, 405.520, 404.830, 404.140, 403.450, 402.760, 402.070, 401.380, 400.690, 
    400.000
]

# Method 1
method1_temps = [
    600.000, 598.462, 596.923, 595.385, 593.846, 592.308, 590.769, 589.231, 587.692, 586.154,
    584.615, 583.077, 581.538, 580.000, 578.461, 576.923, 575.384, 573.846, 572.307, 570.769,
    569.230, 567.692, 566.153, 564.615, 563.077, 561.538, 560.000, 558.461, 556.923, 555.384,
    553.846, 552.307, 550.769, 549.230, 547.692, 546.153, 544.615, 543.077, 541.538, 540.000,
    538.461, 536.923, 535.384, 533.846, 532.307, 530.769, 529.230, 527.692, 526.154, 524.615,
    523.077, 521.538, 520.000, 518.461, 516.923, 515.384, 513.846, 512.308, 510.769, 509.231,
    507.692, 506.154, 504.615, 503.077, 501.538, 500.000, 498.462, 496.923, 495.385, 493.846,
    492.308, 490.769, 489.231, 487.692, 486.154, 484.616, 483.077, 481.539, 480.000, 478.462,
    476.923, 475.385, 473.846, 472.308, 470.770, 469.231, 467.693, 466.154, 464.616, 463.077,
    461.539, 460.000, 458.462, 456.923, 455.385, 453.847, 452.308, 450.770, 449.231, 447.693,
    446.154, 444.616, 443.077, 441.539, 440.000, 438.462, 436.923, 435.385, 433.847, 432.308,
    430.770, 429.231, 427.693, 426.154, 424.616, 423.077, 421.539, 420.000, 418.462, 416.923,
    415.385, 413.846, 412.308, 410.769, 409.231, 407.692, 406.154, 404.615, 403.077, 401.538,
    400.000
]

# Method 2
method2_temps = [
    600.000, 599.540, 599.079, 598.619, 598.159, 597.698, 597.238, 596.778, 596.318, 595.857,
    595.397, 594.258, 592.441, 590.624, 588.807, 586.990, 585.173, 583.356, 581.539, 579.722,
    577.905, 576.088, 574.271, 572.454, 570.637, 568.820, 567.003, 565.186, 563.369, 561.552,
    559.735, 557.918, 556.101, 554.284, 552.467, 550.650, 548.833, 547.016, 545.199, 543.382,
    541.565, 539.748, 537.931, 536.114, 534.297, 532.480, 530.663, 528.846, 527.029, 525.212,
    523.395, 521.578, 519.761, 517.944, 516.127, 514.310, 512.493, 510.676, 508.859, 507.042,
    505.225, 503.408, 501.591, 499.774, 497.957, 496.140, 494.323, 492.506, 490.689, 488.872,
    487.055, 485.238, 483.421, 481.604, 479.787, 477.970, 476.153, 474.336, 472.519, 470.702,
    468.885, 467.068, 465.250, 463.433, 461.616, 459.799, 457.982, 456.165, 454.348, 452.531,
    450.714, 448.897, 447.080, 445.263, 443.446, 441.629, 439.812, 437.995, 436.178, 434.361,
    432.544, 430.726, 428.909, 427.092, 425.275, 423.458, 421.641, 419.824, 418.007, 416.190,
    414.373, 413.119, 412.429, 411.738, 411.048, 410.357, 409.667, 408.976, 408.286, 407.595,
    406.904, 406.213, 405.523, 404.832, 404.142, 403.451, 402.760, 402.070, 401.380, 400.690,
    400.000
]

# Plotting
plt.figure(figsize=(12, 8))

plt.plot(positions, analytical_temps, label='Analytical Solution', color='blue')
plt.plot(positions, method1_temps, label='Method 1', color='red', linestyle='--')
plt.plot(positions, method2_temps, label='Method 2', color='green', linestyle=':')

plt.xlabel('Position (mm)')
plt.ylabel('Temperature (K)')
plt.title('Temperature Distribution Comparison')
plt.legend()
plt.grid(True)

# Save the figure
plt.savefig('4_temperature_distribution_comparison.png', dpi=300)

# Show the plot
plt.show()
