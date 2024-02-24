import { RefObject, useState, useMemo, useEffect } from "react";

export default function useOnScreen(ref: RefObject<HTMLElement>): boolean {
  const [isIntersecting, setIntersecting] = useState<boolean>(false);

  const observer = useMemo(
    () =>
      new IntersectionObserver((entries) => {
        const [entry] = entries;
        // Check if entry is defined before accessing isIntersecting
        if (entry) {
          setIntersecting(entry.isIntersecting);
        }
      }),
    [ref],
  );

  useEffect(() => {
    if (ref.current) {
      observer.observe(ref.current);
    }
    return () => {
      observer.disconnect();
    };
  }, [observer, ref]);

  return isIntersecting;
}
