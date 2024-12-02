using CairoMakie, GLMakie

# num_days = 200
# initial_price = 100.0
# std_dev = 2.0

# price_changes = randn(num_days) * std_dev
# stock_prices = [initial_price]
# for change in price_changes
#     new_price = stock_prices[end] + change
#     push!(stock_prices, new_price)
# end

# fig, ax = lines(1:num_days, stock_prices[1:num_days], color=:red, linestyle=:solid, label="Stock Price")


# function update(frame)
#     lines!(ax, 1:frame, stock_prices[1:frame], color=:red, linestyle=:solid)
#     return ax
# end


# record(fig, "./fig/一维/stock_price_simulation.gif", 1:num_days) do i
#     update(i)
# end

# display(fig)

points = Observable(Point2f[(0, 0)])

fig, ax = scatter(points)
limits!(ax, 0, 30, 0, 30)

frames = 1:30

record(fig, "append_animation.mp4", frames;
        framerate = 30) do frame
    new_point = Point2f(frame, frame)
    points[] = push!(points[], new_point)
end
