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
# Speedup uses the run with Np = 1 and 1 thread at the same N and L.
# --threads-only keeps Np = 1 and varies the OpenMP thread count.

set -o pipefail

root=$(cd "$(dirname "$0")" && pwd)
cd "$root"

# Frontend login shells do not always export Spectrum MPI. The wrapper
# mpixlC and mpisubmit.pl both come from this file on Polus.
if [[ -f /polusfs/setenv/setup.SMPI ]]; then
    # shellcheck disable=SC1091
    source /polusfs/setenv/setup.SMPI
fi

max_jobs=4
do_compile=0
dry_run=0
threads_only=0
threads_wide=0
T_override=""
K=20
poll_sec=15
miss_limit=20
PI=3.141592653589793
results="$root/results"
rows_file=""

usage() {
    cat <<'EOF'
Usage: ./collect_polus.sh [--compile] [--dry-run] [--threads-only] [--threads-wide] [--max-jobs N] [--T value]

  --compile      build ./wave with make (same flags as on Polus)
  --dry-run      print mpisubmit.pl commands and do not submit
  --threads-only one MPI process, OpenMP threads in {1,2,4,8}
                 (Polus allows at most 8 threads on one process)
  --threads-wide 4 MPI processes on 2 nodes, OpenMP threads in {16,32,64}
                 (16 threads use 2 cores each, 32 use 4, 64 use 8)
  --max-jobs N   how many jobs may sit in the queue at once (default 4)
  --T value      one final time for every run; default is the stable step
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
        --threads-only)
            threads_only=1
            shift
            ;;
        --threads-wide)
            threads_wide=1
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
if [[ $threads_only -eq 1 && $threads_wide -eq 1 ]]; then
    die "--threads-only and --threads-wide cannot be used together"
fi
if [[ $threads_only -eq 1 ]]; then
    results="$root/results/threads"
fi
if [[ $threads_wide -eq 1 ]]; then
    results="$root/results/threads-wide"
fi
rows_file="$results/.rows"

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

# Queue "short" rejects a wall-clock limit above 30 minutes.
wtime_for() {
    case $1 in
        128) echo 00:10 ;;
        256) echo 00:20 ;;
        *) echo 00:30 ;;
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
        if [[ $threads_only -eq 1 ]]; then
            # mpisubmit.pl keeps one process only while threads <= 8.
            # A larger -t is rewritten into several MPI processes.
            for threads in 1 2 4 8; do
                queue+=("omp|$n|1|$threads|$ltag|$lval|$tval|$K|$wtime")
            done
        elif [[ $threads_wide -eq 1 ]]; then
            # More than 8 threads per process needs several cores. mpisubmit.pl
            # would also demand 4 nodes, and only 3 compute nodes are up, so
            # these jobs are submitted with bsub: 4 processes, 2 per node.
            for threads in 16 32 64; do
                queue+=("omp|$n|4|$threads|$ltag|$lval|$tval|$K|$wtime")
            done
        else
            for np in 1 4 8 16 32; do
                queue+=("mpi|$n|$np|1|$ltag|$lval|$tval|$K|$wtime")
            done
            for np in 1 2 4 8; do
                queue+=("hybrid|$n|$np|4|$ltag|$lval|$tval|$K|$wtime")
            done
        fi
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
    if [[ $threads_wide -eq 1 ]]; then
        submit_wide=1
        return
    fi
    submit_wide=0
    submit_argv=(
        mpisubmit.pl -p "$np" -t "$threads" -w "$wtime"
        --stdout "$submit_out" --stderr "$submit_err"
        ./wave -- "$n" "$lval" "$lval" "$lval" "$tval" "$k"
    )
}

submit_wide_job() {
    local threads=$1 n=$2 lval=$3 tval=$4 k=$5 wtime=$6
    local slots=$((threads / 2))
    local ptile=$((threads / 4))
    bsub <<EOF
source /polusfs/setenv/setup.SMPI
#BSUB -n ${slots}
#BSUB -R span[ptile=${ptile}]
#BSUB -W ${wtime}
#BSUB -o ${submit_out}
#BSUB -e ${submit_err}
OMP_NUM_THREADS=${threads} mpiexec -npernode 2 ./wave ${n} ${lval} ${lval} ${lval} ${tval} ${k}
EOF
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
        if [[ $submit_wide -eq 1 ]]; then
            echo "bsub -n $(( $(spec_field "$spec" 4) / 2 )) -R span[ptile=$(( $(spec_field "$spec" 4) / 4 ))] -W $(spec_field "$spec" 9) OMP_NUM_THREADS=$(spec_field "$spec" 4) mpiexec -npernode 2 ./wave $(spec_field "$spec" 2) $(spec_field "$spec" 6) $(spec_field "$spec" 6) $(spec_field "$spec" 6) $(spec_field "$spec" 7) $(spec_field "$spec" 8)"
        else
            print_submit
        fi
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
            if (sweep == "wide" && $1 == "omp" && $3 == 4 && $4 == 16 && $12 == "ok" && $9 != "") {
                base[baseline_key($2, $13)] = $9
            } else if ((sweep == "threads" ? $1 == "omp" : $1 == "mpi") && sweep != "wide" && $3 == 1 && $4 == 1 && $12 == "ok" && $9 != "") {
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
            if (sweep == "threads") {
                print "Один процесс. Ускорение S = время(1 нить) / время при тех же N и L." > summary
            } else if (sweep == "wide") {
                print "4 процесса на 2 узлах, по 2 на узел. 16 нитей занимают 2 ядра, 32 — 4, 64 — 8." > summary
                print "Ускорение S = время(16 нитей) / время при тех же N и L." > summary
            } else {
                print "Ускорение S = время(Np=1, 1 нить) / время при тех же N и L." > summary
            }
            print "" > summary

            if (sweep == "threads" || sweep == "wide") {
                ntables = 2
                titles[1] = "OpenMP, L = 1"
                want_mode[1] = "omp"
                want_tag[1] = "L1"
                titles[2] = "OpenMP, L = pi"
                want_mode[2] = "omp"
                want_tag[2] = "Lpi"
            } else {
                ntables = 4
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
            }

            for (table = 1; table <= ntables; table++) {
                print titles[table] > summary
                if (sweep == "wide") {
                    print "Число процессов Np | Число нитей | Число точек N^3 | Время решения | Ускорение S | Погрешность" > summary
                } else if (sweep == "threads") {
                    print "Число нитей | Число точек N^3 | Время решения | Ускорение S | Погрешность" > summary
                } else if (want_mode[table] == "mpi") {
                    print "Число процессов Np | Число точек N^3 | Время решения | Ускорение S | Погрешность" > summary
                } else {
                    print "Число процессов Np | Число нитей | Число точек N^3 | Время решения | Ускорение S | Погрешность" > summary
                }
                nitems = 0
                for (i = 1; i <= nrows; i++) {
                    if (mode[i] != want_mode[table] || ltag[i] != want_tag[table]) {
                        continue
                    }
                    nitems++
                    order[nitems] = i
                }
                for (a = 2; a <= nitems; a++) {
                    held = order[a]
                    b = a - 1
                    while (b >= 1 && (ngrid[order[b]] + 0 > ngrid[held] + 0 || (ngrid[order[b]] + 0 == ngrid[held] + 0 && ((sweep == "threads" || sweep == "wide") ? threads[order[b]] + 0 > threads[held] + 0 : np[order[b]] + 0 > np[held] + 0)))) {
                        order[b + 1] = order[b]
                        b--
                    }
                    order[b + 1] = held
                }
                for (p = 1; p <= nitems; p++) {
                    i = order[p]
                    time_cell = times[i]
                    speed_cell = speed[i]
                    delta_cell = delta[i]
                    if (status[i] != "ok") {
                        time_cell = "fail"
                        speed_cell = "-"
                        delta_cell = "-"
                    }
                    if (sweep == "wide") {
                        print np[i] " | " threads[i] " | " ngrid[i] "^3 | " time_cell " | " speed_cell " | " delta_cell > summary
                    } else if (sweep == "threads") {
                        print threads[i] " | " ngrid[i] "^3 | " time_cell " | " speed_cell " | " delta_cell > summary
                    } else if (want_mode[table] == "mpi") {
                        print np[i] " | " ngrid[i] "^3 | " time_cell " | " speed_cell " | " delta_cell > summary
                    } else {
                        print np[i] " | " threads[i] " | " ngrid[i] "^3 | " time_cell " | " speed_cell " | " delta_cell > summary
                    }
                }
                delete order
                print "" > summary
            }
            close(summary)
        }
    ' results="$results" kshow="$K" toverride="$T_override" sweep="$(
        if [[ $threads_wide -eq 1 ]]; then echo wide
        elif [[ $threads_only -eq 1 ]]; then echo threads
        else echo full
        fi
    )" "$rows_file"
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
    if [[ $submit_wide -eq 1 ]]; then
        lval=$(spec_field "$spec" 6)
        tval=$(spec_field "$spec" 7)
        k=$(spec_field "$spec" 8)
        wtime=$(spec_field "$spec" 9)
        if ! submit_wide_job "$threads" "$n" "$lval" "$tval" "$k" "$wtime" >"$log" 2>&1; then
            cat "$log" >&2
            id=$(sed -n 's/.*[Jj]ob <\([0-9][0-9]*\)>.*/\1/p' "$log" | head -n 1)
            if [[ -n $id ]]; then
                bkill "$id" >/dev/null 2>&1 || true
            fi
            rm -f "$log"
            finish_slot "$spec" "" fail
            return
        fi
        cat "$log"
        id=$(sed -n 's/.*[Jj]ob <\([0-9][0-9]*\)>.*/\1/p' "$log" | head -n 1)
        rm -f "$log"
        if [[ -z $id ]]; then
            echo "bsub did not report a job id for $mode N=$n Np=$np" >&2
            finish_slot "$spec" "" fail
            return
        fi
        slot_id+=("$id")
        slot_spec+=("$spec")
        slot_miss+=(0)
        echo "submit $mode N=$n Np=$np threads=$threads $ltag -> job $id"
        return
    fi
    if ! "${submit_argv[@]}" >"$log" 2>&1; then
        cat "$log" >&2
        id=$(sed -n 's/.*[Jj]ob <\([0-9][0-9]*\)>.*/\1/p' "$log" | head -n 1)
        if [[ -n $id ]]; then
            bkill "$id" >/dev/null 2>&1 || true
        fi
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
