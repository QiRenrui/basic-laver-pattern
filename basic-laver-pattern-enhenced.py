import copy

class BasicLaverPattern:
    def __init__(self, a):
        self.a = copy.deepcopy(a)

    def draw(self, use_fullwidth_space=False):
        base_list = self.a[0][0]
        other_lists = self.a[0][1:]
        max_len = max(max(lst) for lst in other_lists if isinstance(lst, list)) + 1
        result = []

        space_char = '　' if use_fullwidth_space else ' '

        for i, seq in enumerate(other_lists):
            line = [space_char] * max_len
            for num in seq:
                line[num] = 'o'
            if i < len(base_list):
                last_circle_index = max(seq)
                result.append(''.join(line[:last_circle_index + 1]) + f" {base_list[i]}")

        for line in result:
            print(line)

        print("r list:", self.a[1])

    def classify_pattern(self):
        a_0 = self.a[0]
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
        n = len(self.a[0][0])
        self.a[0][0].pop()
        self.a[0].pop()
        if self.a[1]:
            self.a[1].pop()

    def _compute_parent_for_mapped_row(self, row_idx, r_old, u, v, n, b0):
        parent = 0
        if 1 <= row_idx <= len(r_old):
            parent = r_old[row_idx - 1]

        if parent != 0 and u <= parent <= v:
            return parent + (n - u)

        ancestor = parent
        best_candidate = None
        while ancestor != 0:
            if ancestor < b0 and (best_candidate is None or ancestor > best_candidate):
                best_candidate = ancestor
            if 1 <= ancestor <= len(r_old):
                ancestor = r_old[ancestor - 1]
            else:
                break

        if best_candidate is not None:
            return best_candidate

        return 0

    def modify(self):
        a_prime = BasicLaverPattern(copy.deepcopy(self.a))
        a_prime.cut()

        n = len(self.a[0][0])
        l = self.a[0][0][n - 1]
        b = self.a[0][-1][:]
        self.cut()

        u = b[-l - 2]
        v = b[-l - 1] - 1
        self.a[0][0].extend(self.a[0][0][u-1:v])

        r_old = self.a[1][:]
        r_new = r_old[:]
        b0 = b[0]

        for i in range(v - u + 1):
            new_seq = []
            row_idx = u + i
            for elem in self.a[0][row_idx]:
                new_val = self.f(elem, b, u, n, l)
                if new_val == -1:
                    print("Something unpleasant happened. Please contact the author (E-mail: qwerasdfyh@126.com) about the previous pattern so he can improve the rule design.")
                    return a_prime
                new_seq.append(new_val)
            self.a[0].append(sorted(new_seq))
            parent_p = self._compute_parent_for_mapped_row(row_idx, r_old, u, v, n, b0)
            r_new.append(parent_p)

        self.a[1] = r_new

        m = n
        while True:
            if m > len(self.a[0][0]):
                break

            l = self.a[0][0][m - 1]
            e = len(self.a[0][m])

            if e > 2 * l + 1:
                m += 1
                continue

            s = [self.a[0][m][-l - 1]]
            while True:
                s.append(self.a[0][s[-1]][-3])
                if s[-1] <= self.a[0][m][-l - 2]:
                    break

            k = len(s) - 1
            if k == 1:
                m += 1
                continue

            s.pop()
            k -= 1

            for i in range(m, len(self.a[0])):
                self.a[0][i] = [x + k if x >= m + 1 else x for x in self.a[0][i]]

            if e == 2 * l + 1:
                c = self.a[0][m][:-1]
            else:
                c = self.a[0][m][:l - 1] + self.a[0][m][l:-1]

            self.a[0][m].extend(s[1:][::-1] + list(range(m + 1, m + k + 1)))
            self.a[0][m] = sorted(self.a[0][m])
            self.a[0][0][m - 1] += k

            r_before = self.a[1]
            num_rows = len(r_before)

            for j in range(1, num_rows + 1):
                if j == m:
                    continue
                anc = r_before[j - 1]
                while anc != 0:
                    if anc == m:
                        self.a[0][j].extend(s[1:][::-1] + list(range(m + 1, m + k + 1)))
                        self.a[0][j] = sorted(self.a[0][j])
                        self.a[0][0][j - 1] += k
                        break
                    if 1 <= anc <= num_rows:
                        anc = r_before[anc - 1]
                    else:
                        break

            d = []
            for i in range(k):
                d_i = c + s[k - i:] + list(range(m + 1, m + i + 2))
                d.append(sorted(d_i))

            self.a[0][0] = self.a[0][0][:m - 1] + list(range(e - l, e - l + k)) + self.a[0][0][m - 1:]
            self.a[0] = self.a[0][:m] + d + self.a[0][m:]

            old_r = self.a[1]
            R = len(old_r)
            new_r = [0] * (R + k)

            for j in range(1, R + 1):
                parent = old_r[j - 1]
                if j == m:
                    j_new = m
                    parent_new = parent
                else:
                    j_new = j if j < m else j + k
                    if parent != 0 and parent >= m:
                        parent_new = parent + k
                    else:
                        parent_new = parent
                new_r[j_new - 1] = parent_new

            for t in range(1, k + 1):
                idx = m + t
                new_r[idx - 1] = idx - 1

            self.a[1] = new_r

            m += k + 1

        return BasicLaverPattern(self.a)

    def f(self, x, b, u, n, l):
        if x < b[0]:
            return x
        elif x > u:
            return x - u + n
        else:
            for i, val in enumerate(b):
                if x == val:
                    return b[i + l]
            return -1

    def expand(self, n):
        if self.classify_pattern() != "LimitPattern":
            raise ValueError("The pattern must be of LimitPattern type to expand.")
        u = self.a[0][-1][-3]
        n_value = self.a[0][-1][-2]
        d = n_value - u
        b = [u, n_value, n_value + 1]
        self.cut()

        for _ in range(n - 1):
            self.a[0][0].append(1)
            self.a[0].append(b[:])
            self.a[1].append(0)
            self.modify()
            b = [elem + d for elem in b]

        return self


initial_pattern = [
    [
        [1, 1, 2, 2, 2],
        [0, 1, 2],
        [0, 1, 2, 3],
        [0, 1, 2, 3, 4],
        [0, 1, 2, 3, 4, 5],
        [2, 3, 4, 5, 6]
    ],
    [0, 0, 0, 3, 0]
]


def parse_operations(operation_sequence):
    ops = []
    i = 0
    while i < len(operation_sequence):
        if operation_sequence[i] in 'CM':
            ops.append(operation_sequence[i])
            i += 1
        elif operation_sequence[i] == 'E':
            num_str = ''
            i += 1
            while i < len(operation_sequence) and operation_sequence[i].isdigit():
                num_str += operation_sequence[i]
                i += 1
            if not num_str:
                raise ValueError(f"Invalid expand operation: E must be followed by a number.")
            ops.append(f'E{num_str}')
        else:
            raise ValueError(f"Invalid operation character: {operation_sequence[i]}")
    return ops

def reconstruct_pattern_list(operation_sequence):
    pattern_list = [BasicLaverPattern(copy.deepcopy(initial_pattern))]
    try:
        operations = parse_operations(operation_sequence)
    except ValueError as e:
        return None, str(e)

    for op in operations:
        current_pattern = BasicLaverPattern(copy.deepcopy(pattern_list[-1].a))
        if op == 'C':
            if current_pattern.classify_pattern() == "ZeroPattern":
                return None, "Cannot apply cut to a Zero pattern."
            current_pattern.cut()
        elif op == 'M':
            if current_pattern.classify_pattern() != "TransientPattern":
                return None, "Modify operation can only be applied to a Transient pattern."
            current_pattern = current_pattern.modify()
        elif op.startswith('E'):
            try:
                fs_number = int(op[1:])
            except ValueError:
                return None, f"Invalid expand operation: {op}."
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
    def split_operations(sequence):
        ops = []
        i = 0
        while i < len(sequence):
            if sequence[i] in 'MC':
                ops.append(sequence[i])
                i += 1
            elif sequence[i] == 'E':
                num_str = ''
                i += 1
                while i < len(sequence) and sequence[i].isdigit():
                    num_str += sequence[i]
                    i += 1
                ops.append(f'E{num_str}')
        return ops

    s = split_operations(operation_sequence)
    assert len(s) + 1 == len(pattern_list), "Pattern list length doesn't match operation sequence."

    for i in range(len(s)):
        if s[i] == 'E1':
            s[i] = 'C'

    p = [None] * len(pattern_list)
    for i in range(len(s)):
        if s[i] != 'C':
            p[i] = copy.deepcopy(pattern_list[i])
            p[i].cut()
        else:
            p[i] = 0

    i = 0
    while i < len(s):
        if p[i] != 0:
            for j in range(i + 2, len(pattern_list)):
                if p[i].a == pattern_list[j].a:
                    s = s[:i+1] + s[j:]
                    pattern_list = pattern_list[:i+1] + pattern_list[j:]
                    p = p[:i+1] + p[j:]
                    s[i] = 'C'
                    break
        i += 1

    i = 0
    while i < len(s):
        if i >= len(s):
            break
        if s[i].startswith('E'):
            n = int(s[i][1:])
            k = pattern_list[i].a[0][-1][-2] - pattern_list[i].a[0][-1][-3]
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
    pattern_list = [BasicLaverPattern(copy.deepcopy(initial_pattern))]
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

        pattern_type = pattern_list[-1].classify_pattern().replace("Pattern", "").strip()
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
            new_pattern = BasicLaverPattern(copy.deepcopy(pattern_list[-1].a))
            new_pattern.cut()
            pattern_list.append(new_pattern)
            operation_sequence += "C"
            print("Applied cut operation.")

        elif user_input == 'M' and pattern_type == "Transient":
            new_pattern = BasicLaverPattern(copy.deepcopy(pattern_list[-1].a))
            new_pattern = new_pattern.modify()
            pattern_list.append(new_pattern)
            operation_sequence += "M"
            print("Applied modify operation.")

        elif user_input == 'E' and pattern_type == "Limit":
            fs_number = input("Input FS number (positive integer): ").strip()
            if fs_number.isdigit() and int(fs_number) > 0:
                fs_number = int(fs_number)
                new_pattern = BasicLaverPattern(copy.deepcopy(pattern_list[-1].a))
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
                pattern_list = [BasicLaverPattern(copy.deepcopy(initial_pattern))]
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
            input_sequence = input("Input the operation sequence (e.g., MME4C): ").strip().upper()

            if input_sequence == "":
                pattern_list = [BasicLaverPattern(copy.deepcopy(initial_pattern))]
                operation_sequence = ""
            else:
                try:
                    operation_sequence = input_sequence
                    pattern_list, error = reconstruct_pattern_list(operation_sequence)
                    if error:
                        print(f"Error in operation sequence: {error}")
                        pattern_list = [BasicLaverPattern(copy.deepcopy(initial_pattern))]
                        operation_sequence = ""
                except ValueError as e:
                    print(f"Error: {str(e)}")
                    pattern_list = [BasicLaverPattern(copy.deepcopy(initial_pattern))]
                    operation_sequence = ""

        else:
            print("Invalid operation. Please try again.")


main_program()
