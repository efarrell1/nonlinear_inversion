

def debug(*message, gap=True):
    if False:
        print('-----')
        print(*message)
        print('-----')
        if gap:
            print()


def split_string(string, max_length):
    return [string[i:i+max_length] for i in range(0, len(string), max_length)]


def log_parameters(params, gap=True):
    ljust_val = 12
    top_line = [str(x).ljust(ljust_val) for x in params.keys()]
    bottom_line = ["{:.5f}".format(x).ljust(ljust_val) for x in params.values()]

    top_line = '  '.join(top_line)
    bottom_line = '  '.join(bottom_line)

    top_splits = split_string(top_line, (ljust_val + 2)*5)
    bot_splits = split_string(bottom_line, (ljust_val + 2)*5)

    for top, bot in zip(top_splits, bot_splits):
        print(top)
        print(bot)
    
    if gap:
        print()


def log(log_file, *message, gap=True, bold=False):
    if bold:
        message = ["\033[1m" + x + "\033[0m" for x in message]
        print(*message)
    else:
        print(*message)

    if gap:
        print()

    with open(log_file, 'a') as file:
        message = [str(x) for x in message]
        message = ' '.join(message) + '\n'
        if gap:
            message = message + '\n'
        file.write(message)


def log_jacobian(log_file, df_jacob, b, x):
    print("")
    log(log_file, "Jacobian Matrix:", gap=False)
    log(log_file, df_jacob)

    log(log_file, "B vector", ["{:.5f}".format(float(x)) for x in b], gap=False)

    log(log_file, "Solution", x, gap=False)
    # debug("Residuals", residuals, gap=False)
    # debug("Rank", rank)

    log(log_file, "")

def log_initial_status(log_file, nmode_fit, m):
    log(log_file, "Guess for radial order of first mode:", nmode_fit, gap=False)
    log(log_file, "Assumption of l = 1 and m =", m)
    # log(log_file, self.output_dir.split("/")[-1])



def log_separator(log_file):
    log(log_file, "-------------------------------------------------------------------------------------------------------------------------------------")

