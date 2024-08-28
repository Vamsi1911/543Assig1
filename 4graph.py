import matplotlib.pyplot as plt

# Define positions in meters (0 to 0.13 meters with increments)
positions = [i * 0.001 for i in range(131)]

# Temperature data for the first distribution
temp_dist_1 = [
    600.00, 598.50, 597.00, 595.49, 593.98, 592.47, 590.95, 589.43, 587.90, 586.37, 584.84,
    583.30, 581.77, 580.23, 578.69, 577.15, 575.62, 574.08, 572.54, 571.01, 569.47, 567.93,
    566.40, 564.86, 563.32, 561.78, 560.25, 558.71, 557.17, 555.63, 554.10, 552.56, 551.02,
    549.48, 547.95, 546.41, 544.87, 543.33, 541.80, 540.26, 538.72, 537.18, 535.64, 534.11,
    532.57, 531.03, 529.49, 527.95, 526.42, 524.88, 523.34, 521.80, 520.26, 518.73, 517.19,
    515.65, 514.11, 512.57, 511.03, 509.50, 507.96, 506.42, 504.88, 503.34, 501.80, 500.26,
    498.73, 497.19, 495.65, 494.11, 492.57, 491.03, 489.49, 487.95, 486.41, 484.87, 483.34,
    481.80, 480.26, 478.72, 477.18, 475.64, 474.10, 472.56, 471.02, 469.48, 467.94, 466.40,
    464.86, 463.32, 461.78, 460.25, 458.71, 457.17, 455.63, 454.09, 452.55, 451.01, 449.47,
    447.93, 446.39, 444.85, 443.31, 441.77, 440.23, 438.69, 437.15, 435.61, 434.07, 432.53,
    430.99, 429.44, 427.90, 426.36, 424.82, 423.27, 421.73, 420.18, 418.64, 417.09, 415.54,
    413.99, 412.44, 410.89, 409.34, 407.78, 406.23, 404.67, 403.12, 401.56, 400.00
]

# Temperature data for the second distribution
temp_dist_2 = [
    600.00, 598.462, 596.923, 595.385, 593.846, 592.308, 590.769, 589.231, 587.692, 586.154,
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

# Temperature data for the third distribution
temp_dist_3 = [
    600.00, 599.54, 599.079, 598.619, 598.159, 597.698, 597.238, 596.778, 596.318, 595.857,
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
    414.373, 413.119, 412.429, 411.738, 411.048, 410.359, 409.669, 408.980, 408.290
]

# Plotting the temperature distributions
plt.figure(figsize=(12, 8))
plt.plot(positions, temp_dist_1, label='Distribution 1', color='blue')
plt.plot(positions, temp_dist_2, label='Distribution 2', color='green')
plt.plot(positions, temp_dist_3, label='Distribution 3', color='red')
plt.xlabel('Position (m)')
plt.ylabel('Temperature (K)')
plt.title('Comparison of Temperature Distributions')
plt.legend()
plt.grid(True)
plt.show()

