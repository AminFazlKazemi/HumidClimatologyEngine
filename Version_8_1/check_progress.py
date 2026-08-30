import json
import time
from pathlib import Path
import numpy as np

YEAR_DIR = Path(r"C:\c\checkpoints_moisture_v8_0\years")
REFRESH_SECONDS = 30   # هر چند ثانیه یک‌بار به‌روزرسانی شود

def get_progress_snapshot():
    """مجموع واحدهای انجام‌شده و کل واحدها را از فایل‌های JSON موجود می‌خواند."""
    completed = 0
    total = 0
    year_details = []
    for json_file in sorted(YEAR_DIR.glob("year_*.json")):
        try:
            with open(json_file, 'r', encoding='utf-8') as f:
                meta = json.load(f)
            done = int(meta.get('completed_units', 0))
            tot = int(meta.get('total_units', 0))
            year = int(meta.get('year', json_file.stem.split('_')[1]))
            completed += done
            total += tot
            year_details.append((year, done, tot))
        except Exception:
            # فایل خراب یا ناقص را نادیده بگیر
            continue
    return completed, total, year_details

def format_eta(seconds):
    if seconds <= 0:
        return "0h 0m"
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    return f"{hours}h {minutes}m"

def main():
    print("شروع پایش پیشرفت... (برای توقف Ctrl+C را بزنید)")
    prev_completed = None
    prev_time = None

    while True:
        completed, total, details = get_progress_snapshot()
        remaining = total - completed
        pct = 100.0 * completed / total if total > 0 else 0.0

        now = time.time()
        if prev_completed is not None and prev_time is not None:
            dt = now - prev_time
            if dt > 0:
                rate = (completed - prev_completed) / dt   # واحد بر ثانیه
            else:
                rate = 0.0
        else:
            rate = None

        print("\n" + "=" * 60)
        print(f"زمان: {time.strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"پیشرفت کل: {completed:,} از {total:,} واحد ({pct:.2f}%)")
        print(f"باقی‌مانده: {remaining:,} واحد")
        if rate is not None and rate > 0:
            eta_seconds = remaining / rate
            print(f"نرخ فعلی: {rate:.2f} واحد/ثانیه")
            print(f"زمان تقریبی اتمام: {format_eta(eta_seconds)}")
        else:
            print("در حال اندازه‌گیری نرخ... (چند لحظه صبر کنید)")

        # نمایش پیشرفت هر سال (اختیاری، برای سال‌های فعال)
        active_years = [d for d in details if d[1] < d[2]]  # سال‌هایی که هنوز کامل نیستند
        if active_years:
            print("\nسال‌های در حال پردازش:")
            for year, done, tot in active_years[:10]:  # حداکثر ۱۰ سال
                print(f"  {year}: {done}/{tot} ({100.0*done/tot:.1f}%)")

        prev_completed = completed
        prev_time = now

        try:
            time.sleep(REFRESH_SECONDS)
        except KeyboardInterrupt:
            print("\nپایش متوقف شد.")
            break

if __name__ == "__main__":
    main()