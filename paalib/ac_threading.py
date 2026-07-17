def compute_ac_thread_allocation(nthreads, n_amplicons):
    """Divide a thread budget between parallel AC jobs and BFB solver threads.

    BFBArchitect's ILP solvers see diminishing returns beyond roughly two or
    three threads, so multi-amplicon runs favor parallel jobs. A single
    amplicon receives the full thread budget.
    """
    try:
        nthreads = max(1, int(nthreads))
    except (TypeError, ValueError):
        nthreads = 1

    if n_amplicons <= 1:
        return 1, nthreads

    best = None  # (used_threads, jobs, bfb_threads)
    for bfb_threads in (2, 3):
        if bfb_threads > nthreads:
            continue
        jobs = max(1, min(n_amplicons, nthreads // bfb_threads))
        used = jobs * bfb_threads
        # Maximize utilization; on a tie, prefer more parallel jobs.
        if best is None or used > best[0] or (used == best[0] and jobs > best[1]):
            best = (used, jobs, bfb_threads)

    if best is None:  # nthreads == 1
        return 1, 1

    return best[1], best[2]
