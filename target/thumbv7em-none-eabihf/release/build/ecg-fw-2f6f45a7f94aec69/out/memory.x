/* STM32F407VG (Discovery). Под другой камень поменять эти два числа — см. README. */
MEMORY
{
  FLASH : ORIGIN = 0x08000000, LENGTH = 1024K   /* F407VG: 1 МБ flash            */
  RAM   : ORIGIN = 0x20000000, LENGTH = 128K    /* основная SRAM (CCM 64К не трогаем) */
}
