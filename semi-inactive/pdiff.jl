function pdiff(actual::Number, expected::Number)
	return abs((expected - actual)/expected)*100
end

function rdiff(value1::Number, value2::Number, particles::Number)
	return abs((value1 - value2)/particles)*100
end
