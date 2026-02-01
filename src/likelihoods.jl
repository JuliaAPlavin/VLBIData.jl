@stable begin

loglike_parts(model, x::NamedTuple) = loglike_parts(model, x, x.spec)
loglike(model, x::NamedTuple) = loglike(model, x, x.spec)

function loglike_parts(model::F, x::NamedTuple, spec::VisSpec) where {F}
	modvis = visibility(model, spec)
	σ = U.uncertainty(x.value)
	val = U.value(x.value)
	return (modvis, σ, val)
end

function loglike_parts(model::F, x::NamedTuple, spec::VisAmpSpec) where {F}
	modvis = visibility(model, spec)
	μ = abs(modvis)
	σ = U.uncertainty(x.value)
	val = abs(U.value(x.value))
	return (μ, σ, val)
end

function loglike_parts(model::F, x::NamedTuple, spec::ClosurePhaseSpec) where {F}
	modvis = visibility(model, spec)
	μ = angle(modvis)
	aval = angle(x.value)
	σ = U.uncertainty(aval)
	val = U.value(aval)
	return (μ, σ, val)
end

function loglike_parts(model::F, x::NamedTuple, spec::ClosureAmpSpec) where {F}
	modvis = visibility(model, spec)
	μ = log(abs(modvis))
	lval = log(abs(x.value))
	σ = abs(U.uncertainty(lval))
	val = U.value(lval)
	return (μ, σ, val)
end

function loglike(model::F, x::NamedTuple, spec::VisSpec) where {F}
	μ, σ, val = loglike_parts(model, x, spec)
	# equivalent to:
	# logpdf(Normal(real(μ), σ), real(val)) + logpdf(Normal(imag(μ), σ), imag(val))
	# but faster
	two_σ² = 2σ^2
	return -log(π*two_σ²) - abs2(val - μ) / two_σ²
end

function loglike(model::F, x::NamedTuple, spec::VisAmpSpec) where {F}
	μ, σ, val = loglike_parts(model, x, spec)
	# equivalent to logpdf(Rician(μ, σ), val), but without Distributions.jl
	σ² = σ^2
	z = val * μ / σ²
	bx = z < 1e9 ? besselix(0, z) : inv(√(2π * z))
	return log(val / σ² * bx) - (val - μ)^2 / (2σ²)
end

function loglike(model::F, x::NamedTuple, spec::ClosurePhaseSpec) where {F}
	μ, σ, val = loglike_parts(model, x, spec)
	# equivalent to logpdf(VonMises(μ, 1/σ²), val), but without Distributions.jl
	κ = inv(σ^2)
	return κ * (cos(val - μ) - 1) - log(2π * besselix(0, κ))
end

function loglike(model::F, x::NamedTuple, spec::ClosureAmpSpec) where {F}
	μ, σ, val = loglike_parts(model, x, spec)
	!isfinite(μ) && return -Inf
	# equivalent to logpdf(Normal(μ, σ), val) but faster
	two_σ² = 2σ^2
	return -log(π*two_σ²)/2 - (val - μ)^2 / two_σ²
end

end
