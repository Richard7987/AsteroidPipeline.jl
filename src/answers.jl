"""
    answers_dir() -> String

Where the answers of SkyBoT and JPL's Small-Body Identification service are
kept as they came: `\$XDG_CACHE_HOME` (default `~/.cache`)
`/AsteroidPipeline/answers`, one file per query. Every answer is kept;
[`analyze_set`](@ref)`(...; offline=true)` answers from them alone, so an
analysis repeated without the network, or after a service changed its
answers, gives the same result. Deleting it only costs the queries again.
"""
answers_dir() = joinpath(get(ENV, "XDG_CACHE_HOME", joinpath(homedir(), ".cache")), "AsteroidPipeline", "answers")

"""
    OfflineAnswerMissing(service)

Thrown, offline, for a query whose answer was never kept: the set was not
analyzed with the network before, or not with the same frames.
"""
struct OfflineAnswerMissing <: Exception
    service::String
end

Base.showerror(io::IO, e::OfflineAnswerMissing) =
    print(io, "OfflineAnswerMissing: no answer of $(e.service) kept for this query, and working offline ",
          "(analyze the set once with the network)")

const _OFFLINE = Base.ScopedValues.ScopedValue(false)

# `(service, key) -> text` answering in place of the services when set: the tests'
const _ANSWER_RESPONSES = Ref{Any}(nothing)

# The answer to `key` from `service`: asked with `ask()` and kept, if `keep(text)`;
# offline, the answer kept, never asked. `key` must name the query fully and
# nothing that differs between runs.
function _kept_answer(ask, service::AbstractString, key::AbstractString; keep=text -> true)
    dir = joinpath(answers_dir(), service)
    path = joinpath(dir, string(hash(key); base=16) * ".txt")
    if _OFFLINE[]
        isfile(path) || throw(OfflineAnswerMissing(service))
        return read(path, String)
    end
    respond = _ANSWER_RESPONSES[]
    text = respond === nothing ? ask() : respond(service, key)
    if keep(text)
        mkpath(dir)
        tmp = tempname(dir)
        write(tmp, text)
        mv(tmp, path; force=true)
    end
    return text
end

# A query's parameters as one key, in a fixed order
_answer_key(url, query::AbstractDict) = url * "?" * join(["$k=$(query[k])" for k in sort(collect(keys(query)))], "&")
