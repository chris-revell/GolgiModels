
@bind reset Button("Reset")

begin
	reset
	md"Time (t): $(@bind t_end Slider(0.0:0.1:10.0, 10.0, true))"
end

begin
	reset
	md"""
	Initial Rabbit Population:
		$(@bind x_begin Scrubbable(1:5, default = 1)) |
	Initial Fox Population:
		$(@bind y_begin Scrubbable(1:5, default = 1))
	"""
end

begin
	reset
	md"""
	Rabbit Growth:
		$(@bind alpha Scrubbable(0.5:0.01:3.0, default = 2.0, format = ".2f")) |
	Rabbit Decline:
		$(@bind beta Scrubbable(0.5:0.01:3.0, default = 1.2, format = ".2f")) ||
	Fox Growth:
		$(@bind delta Scrubbable(0.5:0.01:3.0, default = 0.9, format = ".2f")) |
	Fox Decline:
		$(@bind gamma Scrubbable(0.5:0.01:3.0, default = 2.9, format = ".2f"))
	"""
end