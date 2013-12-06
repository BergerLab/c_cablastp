char *to_octal_str(int i) {
  char buf[16];
  sprintf(buf, "%o", i);
  return buf;
}
