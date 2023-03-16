import { useState, useEffect } from 'react'

const placeholders = [
    'Paracetamol',
    'Adderall',
    'Alprazolam-Modafinil',
    'Headache',
    'Nausea',
]

const useSearchPlaceholder = (duration) => {
    const [currentPlaceholder, setCurrentPlaceholder] = useState(placeholders[0])
    const [index, setIndex] = useState(0)

    useEffect(() => {
        const interval = setInterval(() => {
            setIndex((prevIndex) => {
                const nextIndex = (prevIndex + 1) % placeholders.length
                setCurrentPlaceholder(placeholders[nextIndex])
                return nextIndex
            })
        }, duration)
        return () => clearInterval(interval)
    }, [duration])

    return currentPlaceholder
}
export default useSearchPlaceholder