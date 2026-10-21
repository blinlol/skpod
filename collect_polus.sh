#!/bin/bash
# Collect wave timings on IBM Polus (frontend polus.hpc.cs.msu.ru).
#
# Build:  http://hpc.cmc.msu.ru/node/241
# Submit: http://hpc.cmc.msu.ru/node/243  (mpisubmit.pl)
# Hybrid: http://hpc.cmc.msu.ru/node/258
# Limits: http://hpc.cmc.msu.ru/node/242  (at most 4 jobs, 200 cores queued)
#
# Default grid is section 6 of the assignment:
#   K = 20
#   L = 1 and L = pi
#   MPI:    Np in {1,4,8,16,32}, 1 OpenMP thread, N in {128,256,512}
#   Hybrid: Np in {1,2,4,8},    4 OpenMP threads, N in {128,256,512}
# Final time T defaults to a stable step, the same for every Np of one (L, N):
#   T = 0.9 * K * L / (N * sqrt(3))
# Speedup uses the MPI run with Np = 1 and 1 thread at the same N and L.

set -o pipefail

root=$(cd "$(dirname "$0")" && pwd)
cd "$root"

max_jobs=4
do_compile=0
dry_run=0
T_override=""
K=20
poll_sec=15
miss_limit=20
PI=3.141592653589793
results="$root/results"
rows_file="$results/.rows"

usage() {
    cat <<'EOF'
Usage: ./collect_polus.sh [--compile] [--dry-run] [--max-jobs N] [--T value]

  --compile    build ./wave with make (same flags as on Polus)
  --dry-run    print mpisubmit.pl commands and do not submit
  --max-jobs N how many jobs may sit in the queue at once (default 4)
  --T value    one final time for every run; default is the stable step
EOF
}

die() {
    echo "$1" >&2
    exit 1
}

positive_number() {
    awk -v x="$1" 'BEGIN {
        ok = x ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ && x + 0 > 0
        exit !ok
    }'
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --compile)
            do_compile=1
            shift
            ;;
        --dry-run)
            dry_run=1
            shift
            ;;
        --max-jobs)
            [[ $# -ge 2 ]] || die "--max-jobs needs a number"
            max_jobs=$2
            shift 2
            ;;
        --max-jobs=*)
            max_jobs=${1#--max-jobs=}
            shift
            ;;
        --T)
            [[ $# -ge 2 ]] || die "--T needs a number"
            T_override=$2
            shift 2
            ;;
        --T=*)
            T_override=${1#--T=}
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            usage >&2
            die "unknown argument: $1"
            ;;
    esac
done

[[ $max_jobs =~ ^[1-9][0-9]*$ ]] || die "--max-jobs must be a positive integer"
if [[ -n $T_override ]]; then
    positive_number "$T_override" || die "--T must be a positive number"
fi

if [[ $do_compile -eq 1 ]]; then
    make || die "make failed"
fi

if [[ $dry_run -eq 0 && ! -x ./wave ]]; then
    die "./wave is not executable; run ./collect_polus.sh --compile or make"
fi

if [[ $dry_run -eq 0 ]]; then
    command -v mpisubmit.pl >/dev/null 2>&1 || die "mpisubmit.pl not found; run this on the Polus frontend"
    command -v bjobs >/dev/null 2>&1 || die "bjobs not found; run this on the Polus frontend"
fi

stable_T() {
    awk -v k="$1" -v l="$2" -v n="$3" 'BEGIN { printf "%.16g", 0.9 * k * l / (n * sqrt(3)) }'
}

wtime_for() {
    case $1 in
        128) echo 00:10 ;;
        256) echo 00:20 ;;
        512) echo 01:00 ;;
        *) echo 01:00 ;;
    esac
}

queue=()
for ltag in L1 Lpi; do
    if [[ $ltag == L1 ]]; then
        lval=1
    else
        lval=$PI
    fi
    for n in 128 256 512; do
        if [[ -n $T_override ]]; then
            tval=$T_override
        else
            tval=$(stable_T "$K" "$lval" "$n")
        fi
        wtime=$(wtime_for "$n")
        for np in 1 4 8 16 32; do
            queue+=("mpi|$n|$np|1|$ltag|$lval|$tval|$K|$wtime")
        done
        for np in 1 2 4 8; do
            queue+=("hybrid|$n|$np|4|$ltag|$lval|$tval|$K|$wtime")
        done
    done
done

spec_field() {
    echo "$1" | cut -d'|' -f"$2"
}

out_path() {
    local mode=$1 n=$2 np=$3 threads=$4 ltag=$5
    echo "$results/${mode}_N${n}_p${np}_t${threads}_${ltag}.out"
}

err_path() {
    local mode=$1 n=$2 np=$3 threads=$4 ltag=$5
    echo "$results/${mode}_N${n}_p${np}_t${threads}_${ltag}.err"
}

fill_submit() {
    local spec=$1
    local mode n np threads ltag lval tval k wtime
    mode=$(spec_field "$spec" 1)
    n=$(spec_field "$spec" 2)
    np=$(spec_field "$spec" 3)
    threads=$(spec_field "$spec" 4)
    ltag=$(spec_field "$spec" 5)
    lval=$(spec_field "$spec" 6)
    tval=$(spec_field "$spec" 7)
    k=$(spec_field "$spec" 8)
    wtime=$(spec_field "$spec" 9)
    submit_out=$(out_path "$mode" "$n" "$np" "$threads" "$ltag")
    submit_err=$(err_path "$mode" "$n" "$np" "$threads" "$ltag")
    submit_argv=(
        mpisubmit.pl -p "$np" -t "$threads" -w "$wtime"
        --stdout "$submit_out" --stderr "$submit_err"
        ./wave -- "$n" "$lval" "$lval" "$lval" "$tval" "$k"
    )
}

print_submit() {
    local arg line
    line=$(printf '%q' "${submit_argv[0]}")
    for arg in "${submit_argv[@]:1}"; do
        line+=$(printf ' %q' "$arg")
    done
    printf '%s\n' "$line" 2>/dev/null || exit 0
}

if [[ $dry_run -eq 1 ]]; then
    for spec in "${queue[@]}"; do
        fill_submit "$spec"
        print_submit
    done
    exit 0
fi

mkdir -p "$results"
: >"$rows_file"

# Bash 4.2 on Polus treats an empty "${array[@]}" under `set -u` as unbound.
# Slots stay indexed from 0; length checks guard every walk.
slot_id=()
slot_spec=()
slot_miss=()

append_row() {
    # mode,N,Np,threads,L,T,K,tau,time_sec,delta,jobid,status,ltag
    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' "$@" >>"$rows_file"
    publish || echo "warning: could not refresh $results/runs.csv" >&2
}

publish() {
    awk -F',' '
        function baseline_key(n, ltag) {
            return n SUBSEP ltag
        }
        function speedup_of(base, time) {
            if (base == "" || time == "" || base + 0 == 0 || time + 0 == 0) {
                return ""
            }
            return sprintf("%.6f", (base + 0) / (time + 0))
        }
        {
            mode[NR] = $1
            ngrid[NR] = $2
            np[NR] = $3
            threads[NR] = $4
            lval[NR] = $5
            tval[NR] = $6
            ksteps[NR] = $7
            tau[NR] = $8
            times[NR] = $9
            delta[NR] = $10
            job[NR] = $11
            status[NR] = $12
            ltag[NR] = $13
            nrows = NR
            if ($1 == "mpi" && $3 == 1 && $4 == 1 && $12 == "ok" && $9 != "") {
                base[baseline_key($2, $13)] = $9
            }
        }
        END {
            csv = results "/runs.csv"
            summary = results "/summary.txt"
            print "mode,N,Np,threads,L,T,K,tau,time_sec,delta,speedup,jobid,status" > csv
            for (i = 1; i <= nrows; i++) {
                s = ""
                if (status[i] == "ok") {
                    s = speedup_of(base[baseline_key(ngrid[i], ltag[i])], times[i])
                }
                print mode[i] "," ngrid[i] "," np[i] "," threads[i] "," lval[i] "," \
                    tval[i] "," ksteps[i] "," tau[i] "," times[i] "," delta[i] "," \
                    s "," job[i] "," status[i] > csv
                speed[i] = s
            }
            close(csv)

            print "Сетка пункта 6. K = " kshow "." > summary
            if (toverride != "") {
                print "T = " toverride " для всех запусков." > summary
            } else {
                print "T = 0.9 * K * L / (N * sqrt(3)), одно значение на пару (L, N)." > summary
            }
            print "Ускорение S = время(Np=1, 1 нить) / время при тех же N и L." > summary
            print "" > summary

            titles[1] = "Таблица 2. MPI, L = 1"
            want_mode[1] = "mpi"
            want_tag[1] = "L1"
            titles[2] = "Таблица 2. MPI, L = pi"
            want_mode[2] = "mpi"
            want_tag[2] = "Lpi"
            titles[3] = "Таблица 3. MPI+OpenMP, L = 1"
            want_mode[3] = "hybrid"
            want_tag[3] = "L1"
            titles[4] = "Таблица 3. MPI+OpenMP, L = pi"
            want_mode[4] = "hybrid"
            want_tag[4] = "Lpi"

            for (table = 1; table <= 4; table++) {
                print titles[table] > summary
                if (want_mode[table] == "mpi") {
                    print "Число процессов Np | Число точек N^3 | Время решения | Ускорение S | Погрешность" > summary
                } else {
                    print "Число процессов Np | Число нитей | Число точек N^3 | Время решения | Ускорение S | Погрешность" > summary
                }
                for (i = 1; i <= nrows; i++) {
                    if (mode[i] != want_mode[table] || ltag[i] != want_tag[table]) {
                        continue
                    }
                    time_cell = times[i]
                    speed_cell = speed[i]
                    delta_cell = delta[i]
                    if (status[i] != "ok") {
                        time_cell = "fail"
                        speed_cell = "-"
                        delta_cell = "-"
                    }
                    if (want_mode[table] == "mpi") {
                        print np[i] " | " ngrid[i] "^3 | " time_cell " | " speed_cell " | " delta_cell > summary
                    } else {
                        print np[i] " | " threads[i] " | " ngrid[i] "^3 | " time_cell " | " speed_cell " | " delta_cell > summary
                    }
                }
                print "" > summary
            }
            close(summary)
        }
    ' results="$results" kshow="$K" toverride="$T_override" "$rows_file"
}

parse_metrics() {
    local file=$1
    local line=""
    parsed_tau=""
    parsed_time=""
    parsed_delta=""
    if [[ ! -f $file ]]; then
        return
    fi
    line=$(grep -E 'time_sec=' "$file" 2>/dev/null | tail -n 1 || true)
    parsed_tau=$(printf '%s\n' "$line" | sed -n 's/.*tau=\([^[:space:]]*\).*/\1/p')
    parsed_time=$(printf '%s\n' "$line" | sed -n 's/.*time_sec=\([^[:space:]]*\).*/\1/p')
    parsed_delta=$(printf '%s\n' "$line" | sed -n 's/.*delta=\([^[:space:]]*\).*/\1/p')
}

record_spec() {
    local spec=$1 jobid=$2 status=$3 tau=$4 time_sec=$5 delta=$6
    local mode n np threads ltag lval tval k
    mode=$(spec_field "$spec" 1)
    n=$(spec_field "$spec" 2)
    np=$(spec_field "$spec" 3)
    threads=$(spec_field "$spec" 4)
    ltag=$(spec_field "$spec" 5)
    lval=$(spec_field "$spec" 6)
    tval=$(spec_field "$spec" 7)
    k=$(spec_field "$spec" 8)
    append_row "$mode" "$n" "$np" "$threads" "$lval" "$tval" "$k" \
        "$tau" "$time_sec" "$delta" "$jobid" "$status" "$ltag"
}

query_job() {
    local id=$1
    local line=""
    line=$(bjobs -a -noheader -o "stat exit_code" "$id" 2>/dev/null | awk 'NF && $1 != "STAT" && $1 != "JOBID" { print; exit }' || true)
    if [[ -n $line ]]; then
        printf '%s\n' "$line"
        return 0
    fi
    line=$(bjobs -a "$id" 2>/dev/null | awk 'NR > 1 && $1 ~ /^[0-9]+$/ { print $3; exit }' || true)
    if [[ -n $line ]]; then
        printf '%s -\n' "$line"
        return 0
    fi
    return 1
}

job_active() {
    case $1 in
        PEND|RUN|PROV|WAIT|PSUSP|USUSP|SSUSP|UNKNOWN) return 0 ;;
        *) return 1 ;;
    esac
}

finish_slot() {
    local spec=$1 jobid=$2 status=$3
    local mode n np threads ltag out
    mode=$(spec_field "$spec" 1)
    n=$(spec_field "$spec" 2)
    np=$(spec_field "$spec" 3)
    threads=$(spec_field "$spec" 4)
    ltag=$(spec_field "$spec" 5)
    out=$(out_path "$mode" "$n" "$np" "$threads" "$ltag")
    parse_metrics "$out"
    if [[ $status == ok ]]; then
        if [[ -z $parsed_time || -z $parsed_delta ]]; then
            status=fail
        fi
    fi
    record_spec "$spec" "$jobid" "$status" "$parsed_tau" "$parsed_time" "$parsed_delta"
    if [[ $status == ok ]]; then
        echo "done job $jobid $mode N=$n Np=$np threads=$threads $ltag time_sec=$parsed_time delta=$parsed_delta"
    else
        echo "fail job ${jobid:-?} $mode N=$n Np=$np threads=$threads $ltag" >&2
    fi
}

submit_one() {
    local spec=$1
    local log id mode n np threads ltag
    mode=$(spec_field "$spec" 1)
    n=$(spec_field "$spec" 2)
    np=$(spec_field "$spec" 3)
    threads=$(spec_field "$spec" 4)
    ltag=$(spec_field "$spec" 5)
    fill_submit "$spec"
    rm -f "$submit_out" "$submit_err"
    log=$(mktemp)
    if ! "${submit_argv[@]}" >"$log" 2>&1; then
        cat "$log" >&2
        rm -f "$log"
        finish_slot "$spec" "" fail
        return
    fi
    cat "$log"
    id=$(sed -n 's/.*[Jj]ob <\([0-9][0-9]*\)>.*/\1/p' "$log" | head -n 1)
    rm -f "$log"
    if [[ -z $id ]]; then
        echo "mpisubmit.pl did not report a job id for $mode N=$n Np=$np" >&2
        finish_slot "$spec" "" fail
        return
    fi
    slot_id+=("$id")
    slot_spec+=("$spec")
    slot_miss+=(0)
    echo "submit $mode N=$n Np=$np threads=$threads $ltag -> job $id"
}

reap() {
    local i nstat ncode spec id miss reply
    local new_id=() new_spec=() new_miss=()
    local count=${#slot_id[@]}
    local mode n np threads ltag out
    for ((i = 0; i < count; i++)); do
        id=${slot_id[$i]}
        spec=${slot_spec[$i]}
        miss=${slot_miss[$i]}
        reply=$(query_job "$id" || true)
        nstat=${reply%% *}
        ncode=${reply##* }
        if [[ -z $reply ]]; then
            mode=$(spec_field "$spec" 1)
            n=$(spec_field "$spec" 2)
            np=$(spec_field "$spec" 3)
            threads=$(spec_field "$spec" 4)
            ltag=$(spec_field "$spec" 5)
            out=$(out_path "$mode" "$n" "$np" "$threads" "$ltag")
            parse_metrics "$out"
            if [[ -n $parsed_time && -n $parsed_delta ]]; then
                finish_slot "$spec" "$id" ok
                continue
            fi
            miss=$((miss + 1))
            if ((miss >= miss_limit)); then
                finish_slot "$spec" "$id" fail
                continue
            fi
            new_id+=("$id")
            new_spec+=("$spec")
            new_miss+=("$miss")
            continue
        fi
        if job_active "$nstat"; then
            new_id+=("$id")
            new_spec+=("$spec")
            new_miss+=(0)
            continue
        fi
        if [[ $nstat == DONE && ( $ncode == 0 || $ncode == "-" || $ncode == DONE ) ]]; then
            mode=$(spec_field "$spec" 1)
            n=$(spec_field "$spec" 2)
            np=$(spec_field "$spec" 3)
            threads=$(spec_field "$spec" 4)
            ltag=$(spec_field "$spec" 5)
            out=$(out_path "$mode" "$n" "$np" "$threads" "$ltag")
            parse_metrics "$out"
            if [[ -n $parsed_time && -n $parsed_delta ]]; then
                finish_slot "$spec" "$id" ok
                continue
            fi
            miss=$((miss + 1))
            if ((miss >= miss_limit)); then
                finish_slot "$spec" "$id" fail
                continue
            fi
            new_id+=("$id")
            new_spec+=("$spec")
            new_miss+=("$miss")
            continue
        fi
        finish_slot "$spec" "$id" fail
    done
    if ((${#new_id[@]} == 0)); then
        slot_id=()
        slot_spec=()
        slot_miss=()
    else
        slot_id=("${new_id[@]}")
        slot_spec=("${new_spec[@]}")
        slot_miss=("${new_miss[@]}")
    fi
}

idx=0
total=${#queue[@]}
while ((idx < total)) || ((${#slot_id[@]} > 0)); do
    reap
    while ((idx < total)) && ((${#slot_id[@]} < max_jobs)); do
        submit_one "${queue[$idx]}"
        idx=$((idx + 1))
    done
    if ((${#slot_id[@]} > 0)); then
        sleep "$poll_sec"
    fi
done

echo "results: $results/runs.csv"
echo "summary: $results/summary.txt"
