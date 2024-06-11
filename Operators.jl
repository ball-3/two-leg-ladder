module Operators
#TODO souper douper unfiniished

const AS_MPO = 1
const AS_OM = 2

const HAM = HAMILTONIAN = 1
const SPIN_P_OP = SPIN_PLUS_OPERATOR = 2
const SPIN_M_OP = SPIN_MINUS_OPERATOR = 3
const SPIN_Z_OP = SPIN_Z_OPERATOR = 4
const SPIN_X_OP = SPIN_X_OPERATOR = 5
const SPIN_Y_OP = SPIN_Y_OPERATOR = 6

formatError = ArgumentError("format specified was not in the list of known formats")
operatorError = ArgumentError("operator specified was not in the list of known operators")

#taking format to be AS_MPO or AS_OM
#taking operator to be one of the defined constants 1-6
#indices should be a sites array if MPO, and number 
function getOperator(format,operator,indices,indiciesAppliedTo)
	if operator == 1
		if format == 1
		elseif format == 2
		else
			return formatError
		end
	elseif operator == 2
		if format == 1
		elseif format == 2
		else
			return formatError
		end
	elseif operator == 3
		if format == 1
		elseif format == 2
		else
			return formatError
		end
	elseif operator == 4
		if format == 1
		elseif format == 2
		else
			return formatError
		end
	elseif operator == 5
		if format == 1
		elseif format == 2
		else
			return formatError
		end
	elseif operator == 6
	else
		return operatorError
	end
end

function spinOM(particles, locations...)
	
	len = length(locations)

	if len == 0
		dim = 2^particles
		return zeros(Complex{Float64},dim,dim)
	elseif len == 1
		#return matrices[1]
	end

end

function spinMPO()
end

end
