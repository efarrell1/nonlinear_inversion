

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


def log_jacobian(df_jacob, b, x):
    print("")
    log(self.log_file, "Jacobian Matrix:", gap=False)
    log(self.log_file, self.df_jacob)

    log(self.log_file, "B vector", ["{:.5f}".format(float(x)) for x in self.b], gap=False)

    log(self.log_file, "Solution", self.x, gap=False)
    # debug("Residuals", residuals, gap=False)
    # debug("Rank", rank)

    log(self.log_file, "")




def log_separator(log_file):
    log(log_file, "-------------------------------------------------------------------------------------------------------------------------------------")

