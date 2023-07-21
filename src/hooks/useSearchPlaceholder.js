import { useState, useEffect } from 'react';
import { searchTypes } from "../components/OptionModal/OptionModal";

const placeholders = {
    [searchTypes[0].value]: ['Acetaminophen', 'Alprazolam', 'Oxycodone', 'Modafinil', 'Ibuprofen', 'Diazepam',
        'Hydrocodone', 'Tramadol', 'Codeine', 'Gabapentin', 'Meloxicam', 'Cyclobenzaprine', 'Naproxen',
        'Methocarbamol', 'Prednisone', 'Citalopram', 'Amitriptyline', 'Trazodone', 'Lisinopril', 'Atorvastatin',
        'Metformin', 'Amlodipine', 'Omeprazole', 'Metoprolol', 'Simvastatin', 'Losartan', 'Azithromycin',
        'Hydrochlorothiazide', 'Amoxicillin', 'Albuterol', 'Levothyroxine', 'Furosemide', 'Fluticasone',
        'Montelukast', 'Fluoxetine', 'Escitalopram', 'Sertraline', 'Bupropion', 'Cyclobenzaprine', 'Trazodone',
        'Lisinopril', 'Atorvastatin', 'Metformin', 'Amlodipine', 'Omeprazole', 'Metoprolol', 'Simvastatin', 'Losartan',
        'Azithromycin', 'Hydrochlorothiazide', 'Amoxicillin', 'Albuterol', 'Levothyroxine', 'Furosemide',
        'Fluticasone', 'Montelukast', 'Fluoxetine', 'Escitalopram', 'Sertraline', 'Bupropion'],
    [searchTypes[1].value]: ['Xanax', 'Percocet', 'Adderall', 'Valium', 'Vicodin', 'Ambien', 'Klonopin', 'Oxycontin',
        'Concerta', 'Ritalin', 'Zoloft', 'Ativan', 'Lyrica', 'Lunesta', 'Lexapro', 'Prozac'],
    [searchTypes[2].value]: ['Headache', 'Nausea', 'Fever'],
};

const useSearchPlaceholder = (duration, searchType) => {
    const [currentPlaceholder, setCurrentPlaceholder] = useState('');
    const [currentIndex, setCurrentIndex] = useState(null);

    useEffect(() => {
        const setRandomPlaceholder = () => {
            let randomIndex;
            do {
                randomIndex = Math.floor(Math.random() * placeholders[searchType].length);
            } while (randomIndex === currentIndex);

            setCurrentPlaceholder(placeholders[searchType][randomIndex]);
            setCurrentIndex(randomIndex);
        };

        setRandomPlaceholder(); // set a random placeholder immediately

        const interval = setInterval(setRandomPlaceholder, duration); // then update it every `duration` milliseconds

        return () => clearInterval(interval);
    }, [duration, searchType]) // eslint-disable-line react-hooks/exhaustive-deps

    return currentPlaceholder;
};

export default useSearchPlaceholder;
