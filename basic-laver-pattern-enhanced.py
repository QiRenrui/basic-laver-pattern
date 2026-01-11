'''
This version of blp do not guarantee that there's no infinite descending chain. It's not hard to guarantee it by verify some intermediate informations during the modify operation, but for standard patterns it will never be used, it will only cause more time in computation.
'''
import bisect

def _clone_rows(rows):
    return [row[:] for row in rows]

def _clone_mask(mask):
    return [set(s) for s in mask]

def _find_index(sorted_row, val):
    i = bisect.bisect_left(sorted_row, val)
    if i < len(sorted_row) and sorted_row[i] == val:
        return i
    return None

def _shift_sorted_row_inplace(sorted_row, threshold, delta):
    if delta == 0:
        return
    i = bisect.bisect_left(sorted_row, threshold)
    for j in range(i, len(sorted_row)):
        sorted_row[j] += delta

def _shift_mark_set(mark_set, threshold, delta):
    if delta == 0 or not mark_set:
        return mark_set
    new = set()
    for x in mark_set:
        new.add(x + delta if x >= threshold else x)
    return new


class BasicLaverPattern:
    def __init__(self, rows, mask=None):
        self.rows = _clone_rows(rows)
        if mask is None:
            self.mask = [set() for _ in self.rows]
        else:
            self.mask = _clone_mask(mask)
        self.mask[0] = set()

    def clone(self):
        return BasicLaverPattern(self.rows, self.mask)

    def draw(self, use_fullwidth_space=False):
        base_list = self.rows[0]
        other_lists = self.rows[1:]
        if not other_lists:
            print("(empty)")
            return
        max_len = max(seq[-1] for seq in other_lists) + 1
        space_char = '　' if use_fullwidth_space else ' '
        result = []
        for i, seq in enumerate(other_lists, start=1):
            line = [space_char] * max_len
            mset = self.mask[i]
            for num in seq:
                line[num] = 'a' if num in mset else 'o'
            if i <= len(base_list):
                last_circle_index = seq[-1]
                result.append(''.join(line[:last_circle_index + 1]) + f" {base_list[i-1]}")
        for line in result:
            print(line)

    def classify_pattern(self):
        a_0 = self.rows
        if len(a_0[0]) == 2:
            return "ZeroPattern"
        elif len(a_0[-1]) == 5 and a_0[-1][2] == 2:
            return "SuccessorPattern"
        elif len(a_0[-1]) == 6 and a_0[0][-1] == 3 and a_0[-1][1] == 1:
            return "LimitPattern"
        else:
            return "TransientPattern"

    def cut(self):
        if self.classify_pattern() == "ZeroPattern":
            print("This is a ZeroPattern. The cut operation will not be applied.")
            return
        self.rows[0].pop()
        self.rows.pop()
        self.mask.pop()

    def _transmission_penultimate_and_terminal(self, row_idx, n_value):
        rows = self.rows
        base = rows[0]
        row = rows[row_idx]
        l_m = base[row_idx - 1]
        k = _find_index(row, n_value)
        if k is None or k - l_m < 0:
            raise RuntimeError("Unpleasant.")
        threshold = row[k - l_m]
        cur = n_value
        prev = None
        visited = {cur}
        while True:
            l_s = base[cur - 1]
            nxt = rows[cur][-l_s - 2]
            prev, cur = cur, nxt
            if cur <= threshold:
                return prev, cur
            if cur in visited:
                raise RuntimeError("Unpleasant.")
            visited.add(cur)

    def _first_not_copied_in_transmission(self, orig_rows, orig_base, copied_set, row_idx, n_value):
        row = orig_rows[row_idx]
        l_m = orig_base[row_idx - 1]
        k = _find_index(row, n_value)
        if k is None or k - l_m < 0:
            raise RuntimeError("Unpleasant.")
        threshold = row[k - l_m]
        cur = n_value
        seq = [cur]
        visited = {cur}
        while True:
            l_s = orig_base[cur - 1]
            nxt = orig_rows[cur][-l_s - 2]
            seq.append(nxt)
            cur = nxt
            if cur <= threshold:
                break
            if cur in visited:
                raise RuntimeError("Unpleasant.")
            visited.add(cur)
        t = seq[-2]
        terminal = seq[-1]
        tprime = None
        for x in seq:
            if x not in copied_set:
                tprime = x
                break
        return tprime, t, terminal

    def _slice_right_block(self, row_idx, anchor, q):
        row = self.rows[row_idx]
        pos = _find_index(row, anchor)
        if pos is None:
            raise RuntimeError("Unpleasant.")
        block = row[pos + 1: pos + 1 + q]
        if len(block) != q:
            raise RuntimeError("Unpleasant.")
        return block

    def _mark_completion_for_row(self, r, meta, native_done):
        base = self.rows[0]
        row0 = self.rows[r]
        initial_marks = [x for x in row0 if x in self.mask[r]]
        before = set(row0)
        added_total = 0
        for n in initial_marks:
            if n not in self.mask[r]:
                continue
            if _find_index(self.rows[r], n) is None:
                continue
            if n <= 0 or n >= len(meta):
                continue
            info = meta[n]
            if not info or not info.get("native_generated", False):
                continue
            t, n_terminal = self._transmission_penultimate_and_terminal(r, n)
            q = native_done.get(t, 0)
            if q <= 0:
                continue
            target_row = t + q
            left_block = self._slice_right_block(target_row, n_terminal, q)
            right_block = list(range(n + 1, n + q + 1))
            new_vals = set(left_block) | set(right_block)
            truly_new = new_vals - before
            if truly_new:
                added_total += len(truly_new)
                before |= truly_new
            row_set = set(self.rows[r])
            row_set.update(new_vals)
            self.rows[r] = sorted(row_set)
            self.mask[r].difference_update(left_block)
            self.mask[r].update(right_block)
        if added_total > 0:
            base[r - 1] += (added_total // 2)

    def _shift_values_ge(self, start_row_idx, threshold, delta):
        for i in range(start_row_idx, len(self.rows)):
            _shift_sorted_row_inplace(self.rows[i], threshold, delta)
            self.mask[i] = _shift_mark_set(self.mask[i], threshold, delta)

    def _native_completion_step(self, m, meta):
        rows = self.rows
        base = rows[0]
        l = base[m - 1]
        e = len(rows[m])
        if e > 2 * l + 1:
            return False, 0
        s = [rows[m][-l - 1]]
        while True:
            s.append(rows[s[-1]][-3])
            if s[-1] <= rows[m][-l - 2]:
                break
        k = len(s) - 1
        if k == 1:
            return False, 0
        s.pop()
        q = k - 1
        if q <= 0:
            return False, 0
        marks_m_orig = set(self.mask[m])
        self._shift_values_ge(m, m + 1, q)
        if e == 2 * l + 1:
            c = rows[m][:-1]
        else:
            c = rows[m][:l - 1] + rows[m][l:-1]
        ext = s[1:][::-1] + list(range(m + 1, m + q + 1))
        rows[m].extend(ext)
        rows[m].sort()
        base[m - 1] += q
        d = []
        for i in range(q):
            d_i = c + s[q - i:] + list(range(m + 1, m + i + 2))
            d.append(sorted(d_i))
        base[:] = base[:m - 1] + list(range(e - l, e - l + q)) + base[m - 1:]
        rows[:] = rows[:m] + d + rows[m:]
        self.mask[:] = self.mask[:m] + [set() for _ in range(q)] + self.mask[m:]
        meta_insert = [{"native_generated": True, "native_q": q} for _ in range(q)]
        meta[:] = meta[:m] + meta_insert + meta[m:]
        marks_to_propagate = {x + q if x >= m + 1 else x for x in marks_m_orig}
        for row_idx in range(m, m + q + 1):
            self.mask[row_idx].update(marks_to_propagate)
        for j in range(1, q + 1):
            self.mask[m + j].update(range(m, m + j))
        self.mask[m + q].discard(m + 1 + q)
        return True, q

    def modify(self, copy_only=False):
        a_prime = self.clone()
        a_prime.cut()
        orig_rows = _clone_rows(self.rows)
        orig_mask = _clone_mask(self.mask)
        orig_base = orig_rows[0][:]
        orig_last_mask = set(orig_mask[-1])
        rows = self.rows
        base0 = rows[0]
        n_before_cut = len(base0)
        l_last = base0[n_before_cut - 1]
        b = rows[-1][:]
        b0 = b[0]
        p_leftmost = b[0]
        self.cut()
        rows = self.rows
        base0 = rows[0]
        u = b[-l_last - 2]
        v = b[-l_last - 1] - 1
        base0.extend(base0[u - 1: v])
        b_map = {}
        limit = len(b) - l_last
        for i in range(limit):
            key = b[i]
            if key not in b_map:
                b_map[key] = b[i + l_last]
        def map_elem(x):
            if x < b0:
                return x
            if x > u:
                return x - u + n_before_cut
            return b_map.get(x, -1)
        copied_set = set(range(u, v + 1))
        for i in range(v - u + 1):
            row_idx = u + i
            src_row = rows[row_idx]
            new_seq = []
            for elem in src_row:
                new_val = map_elem(elem)
                if new_val == -1:
                    print("Something unpleasant happened. Please contact the author (E-mail: qwerasdfyh@126.com) about the previous pattern so he can improve the rule design.")
                    return a_prime
                new_seq.append(new_val)
            new_seq.sort()
            rows.append(new_seq)
            new_marks = set()
            src_marks = orig_mask[row_idx]
            if src_marks:
                l_m = orig_base[row_idx - 1]
                for marked_val in src_marks:
                    if _find_index(orig_rows[row_idx], marked_val) is None:
                        continue
                    tprime, t, _terminal = self._first_not_copied_in_transmission(
                        orig_rows, orig_base, copied_set, row_idx, marked_val
                    )
                    keep = False
                    if t in copied_set:
                        keep = True
                    else:
                        if tprime is not None and tprime < p_leftmost:
                            keep = True
                        elif tprime is not None:
                            u_img = b_map.get(tprime, None)
                            if u_img is not None and u_img in orig_last_mask:
                                pos_u = _find_index(new_seq, u_img)
                                if pos_u is not None:
                                    idx_check = pos_u - l_m + 1
                                    if idx_check >= 0 and new_seq[idx_check] <= p_leftmost:
                                        keep = True
                    if keep:
                        new_marks.add(map_elem(marked_val))
            self.mask.append(new_marks)
        if copy_only:
            return self.clone()
        meta = [None] * len(self.rows)
        native_done = {}
        m = n_before_cut
        while True:
            base0 = self.rows[0]
            if m > len(base0):
                break
            self._mark_completion_for_row(m, meta, native_done)
            did, q = self._native_completion_step(m, meta)
            if did:
                native_done[m] = q
                m += q + 1
            else:
                m += 1
        return self.clone()

    def expand(self, n):
        if self.classify_pattern() != "LimitPattern":
            raise ValueError("The pattern must be of LimitPattern type to expand.")
        u = self.rows[-1][-3]
        n_value = self.rows[-1][-2]
        d = n_value - u
        b = [u, n_value, n_value + 1]
        self.cut()
        for _ in range(n - 1):
            self.rows[0].append(1)
            self.rows.append(b[:])
            self.mask.append(set())
            self.modify(copy_only=True)
            b = [elem + d for elem in b]
        return self


initial_rows = [
    [1, 1, 2, 2, 2],
    [0, 1, 2],
    [0, 1, 2, 3],
    [0, 1, 2, 3, 4],
    [0, 1, 2, 3, 4, 5],
    [2, 3, 4, 5, 6]
]
initial_mask = [set() for _ in initial_rows]
initial_mask[4] = {3}


def parse_operations(operation_sequence):
    ops = []
    i = 0
    while i < len(operation_sequence):
        c = operation_sequence[i]
        if c in 'CM':
            ops.append(c)
            i += 1
        elif c == 'E':
            num_str = ''
            i += 1
            while i < len(operation_sequence) and operation_sequence[i].isdigit():
                num_str += operation_sequence[i]
                i += 1
            if not num_str:
                raise ValueError("Invalid expand operation: E must be followed by a number.")
            ops.append(f'E{num_str}')
        else:
            raise ValueError(f"Invalid operation character: {operation_sequence[i]}")
    return ops

def reconstruct_pattern_list(operation_sequence):
    pattern_list = [BasicLaverPattern(initial_rows, initial_mask)]
    try:
        operations = parse_operations(operation_sequence)
    except ValueError as e:
        return None, str(e)
    for op in operations:
        current_pattern = pattern_list[-1].clone()
        if op == 'C':
            if current_pattern.classify_pattern() == "ZeroPattern":
                return None, "Cannot apply cut to a Zero pattern."
            current_pattern.cut()
        elif op == 'M':
            if current_pattern.classify_pattern() != "TransientPattern":
                return None, "Modify operation can only be applied to a Transient pattern."
            current_pattern = current_pattern.modify()
        elif op.startswith('E'):
            fs_number = int(op[1:])
            if fs_number < 1:
                return None, "Expand number must be a positive integer."
            if current_pattern.classify_pattern() != "LimitPattern":
                return None, "Expand operation can only be applied to a Limit pattern."
            current_pattern.expand(fs_number)
        else:
            return None, f"Invalid operation: {op}."
        pattern_list.append(current_pattern)
    return pattern_list, None

def optimize_simplify(operation_sequence, pattern_list):
    s = parse_operations(operation_sequence)
    assert len(s) + 1 == len(pattern_list)
    for i in range(len(s)):
        if s[i] == 'E1':
            s[i] = 'C'
    p = [None] * len(pattern_list)
    for i in range(len(s)):
        if s[i] != 'C':
            p[i] = pattern_list[i].clone()
            p[i].cut()
        else:
            p[i] = 0
    i = 0
    while i < len(s):
        if p[i] != 0:
            for j in range(i + 2, len(pattern_list)):
                if p[i].rows == pattern_list[j].rows and p[i].mask == pattern_list[j].mask:
                    s = s[:i+1] + s[j:]
                    pattern_list = pattern_list[:i+1] + pattern_list[j:]
                    p = p[:i+1] + p[j:]
                    s[i] = 'C'
                    break
        i += 1
    i = 0
    while i < len(s):
        if s[i].startswith('E'):
            n = int(s[i][1:])
            k = pattern_list[i].rows[-1][-2] - pattern_list[i].rows[-1][-3]
            if i + k < len(s) and all(x == 'C' for x in s[i+1:i+k+1]):
                s[i] = f'E{n-1}'
                s = s[:i+1] + s[i+k+1:]
                pattern_list = pattern_list[:i+1] + pattern_list[i+k+1:]
            else:
                i += 1
        else:
            i += 1
    operation_list = ''.join(s)
    return operation_list, pattern_list


def main_program():
    pattern_list = [BasicLaverPattern(initial_rows, initial_mask)]
    operation_sequence = ""
    SpaceWidth = False
    while True:
        if len(pattern_list) == 1:
            print("""
This is an ordinal notation analysis tool for exploring the lower bounds of numbers in the Laver table.
Author's E-mail address: qwerasdfyh@126.com.

Each pattern represents an ordinal. We can define a function, f(p, n), where p is a pattern and n is a positive integer:

- If p is a zero pattern, then f(p, n) = n^n.
- If p is a transient pattern that can be modified to q, then f(p, n) = f(q, n).
- If p is a successor pattern and will be cut to q, then f(p, n) is the result of iterating n using f(q, _ ) for n times.
- If p is a limit pattern, then f(p, n) = f(p[n], n).

For the initial pattern p_init, the author has proven f(p_init, n) is well-defined, and γ_F(4) > j₍₁₀₎(κ₃) ≥ γ_f(p_init, 2), where F is the Laver function, providing a lower bound for F(4), that is much stronger than Dougherty's 1992 bound.
The author also proved that F(n+3) > f(p_init, 2^n), for any positive integer n. There should be a lot of space for this bound to improve.
            """)
        print("\nCurrent pattern:")
        pattern_list[-1].draw(use_fullwidth_space=SpaceWidth)
        print(f"Operation sequence: {operation_sequence}")
        pclass = pattern_list[-1].classify_pattern()
        pattern_type = pclass.replace("Pattern", "").strip()
        message = f"This is a {pattern_type} pattern."
        if pattern_type != "Zero":
            message += " C: Cut."
        if pattern_type == "Transient":
            message += " M: Modify."
        elif pattern_type == "Limit":
            message += " E: Expand."
        if len(pattern_list) > 1:
            message += " U: Undo."
        message += " S: Simplify."
        message += " I: Input operations."
        message += " W: SpaceWidth toggle."
        print(message)
        user_input = input("Enter your operation: ").strip().upper()
        if user_input == 'W':
            SpaceWidth = not SpaceWidth
            print(f"Space width toggled to {'Full-width' if SpaceWidth else 'Half-width'}.")
        elif user_input == 'C' and pattern_type != "Zero":
            new_pattern = pattern_list[-1].clone()
            new_pattern.cut()
            pattern_list.append(new_pattern)
            operation_sequence += "C"
            print("Applied cut operation.")
        elif user_input == 'M' and pattern_type == "Transient":
            new_pattern = pattern_list[-1].clone()
            new_pattern = new_pattern.modify()
            pattern_list.append(new_pattern)
            operation_sequence += "M"
            print("Applied modify operation.")
        elif user_input == 'E' and pattern_type == "Limit":
            fs_number = input("Input FS number (positive integer): ").strip()
            if fs_number.isdigit() and int(fs_number) > 0:
                fs_number = int(fs_number)
                new_pattern = pattern_list[-1].clone()
                new_pattern.expand(fs_number)
                pattern_list.append(new_pattern)
                operation_sequence += f"E{fs_number}"
                print(f"Expanded the pattern {fs_number} times.")
            else:
                print("Invalid FS number. Please enter a positive integer.")
        elif user_input == 'U' and len(pattern_list) > 1:
            operations = parse_operations(operation_sequence)
            operations = operations[:-1]
            operation_sequence = ''.join(operations)
            pattern_list, error = reconstruct_pattern_list(operation_sequence)
            if error:
                print(f"Error in operation sequence: {error}")
                pattern_list = [BasicLaverPattern(initial_rows, initial_mask)]
                operation_sequence = ""
            print("Undo the last operation.")
        elif user_input == 'S':
            optimized_sequence, optimized_pattern_list = optimize_simplify(operation_sequence, pattern_list)
            if optimized_sequence != operation_sequence:
                print(f"Simplified operation sequence: {optimized_sequence}")
                operation_sequence = optimized_sequence
                pattern_list = optimized_pattern_list
            else:
                print("No further simplifications possible.")
        elif user_input == 'I':
            input_sequence = input("Input the operation sequence (e.g., MCCE2MMCCE2MCMCC): ").strip().upper()
            if input_sequence == "":
                pattern_list = [BasicLaverPattern(initial_rows, initial_mask)]
                operation_sequence = ""
            else:
                operation_sequence = input_sequence
                pattern_list, error = reconstruct_pattern_list(operation_sequence)
                if error:
                    print(f"Error in operation sequence: {error}")
                    pattern_list = [BasicLaverPattern(initial_rows, initial_mask)]
                    operation_sequence = ""
        else:
            print("Invalid operation. Please try again.")


main_program()
